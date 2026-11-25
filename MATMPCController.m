classdef MATMPCController < handle
    %MATMPCCONTROLLER Run one instance of a built MATMPC model.
    % C = MATMPCController(SETTINGS, OPTIONS) initializes a built example.
    % C = MATMPCController(SETTINGS, OPTIONS, INPUT) supplies its initial data.
    % [U, SOLUTION] = C.step(X0, UPDATES) applies numeric runtime data and solves.
    % Updates to references, weights, parameters and bounds persist until reset
    % or another update for the same field.
    % C.snapshot(true) includes copies of the last QP workspace arrays.
    % C.reset() starts again from the initial data and releases solver state.
    % Multiple instances share one active model and MEX set. Delete them
    % before changing the model, MATLAB path, or loaded solver MEX files.

    properties (Access = private)
        Settings
        Options
        Initial
        Input
        Memory
        Info
        Lease
        Sample = 0
    end

    methods
        function obj = MATMPCController(settings, opt, input)
            if nargin < 2
                opt = struct();
            end
            project = fileparts(mfilename('fullpath'));
            obj.Options = controllerOptions(opt);
            settings = controllerSettings(project, settings, obj.Options);
            obj.Settings = settings;

            if nargin < 3
                if obj.Options.nonuniform_grid
                    [input, data] = InitData_ngrid(settings);
                else
                    [input, data] = InitData(settings);
                end
                reference = data.REF(1, :).';
                input.y = repmat(reference, 1, settings.N);
                input.yN = reference(1:settings.nyN);
            end
            obj.Initial = validateInput(settings, input);
            obj.reset();
            obj.Lease = onCleanup(@() controllerCount(-1));
            controllerCount(1);
        end

        function [u, solution] = step(obj, measurement, updates)
            validateattributes(measurement, {'double'}, ...
                {'real', 'finite', 'nonsparse', 'vector', 'numel', obj.Settings.nx}, ...
                mfilename, 'measurement');
            input = obj.Input;
            input.x0 = measurement(:);
            if nargin > 2
                input = applyUpdates(obj.Settings, input, updates);
            end
            if isempty(obj.Memory)
                obj.Memory = obj.initializeMemory();
            end
            mem = obj.Memory;
            mem.iter = obj.Sample + 1;
            try
                [output, mem] = mpc_nmpcsolver(input, obj.Settings, mem, obj.Options);
            catch problem
                % Legacy kernels may have overwritten some of the workspace.
                % Discard it and keep the last committed MATLAB iterate.
                try
                    obj.releaseMemory();
                catch cleanupProblem
                    problem = addCause(problem, cleanupProblem);
                end
                rethrow(problem);
            end

            obj.Memory = mem;
            obj.Sample = mem.iter;
            obj.Info = output.info;
            obj.Input = copyFields(input, output, iterateNames());
            u = output.u(:, 1);
            if nargout > 1
                solution = copyNumericFields(output, iterateNames());
                solution.info = output.info;
            end
        end

        function state = snapshot(obj, includeQP)
            if nargin < 2
                includeQP = false;
            end
            validateattributes(includeQP, {'logical'}, {'scalar'}, mfilename, 'includeQP');
            state.settings = obj.Settings;
            state.options = obj.Options;
            state.sample = obj.Sample;
            state.input = copyNumericFields(obj.Input, fieldnames(obj.Input));
            state.info = obj.Info;

            if includeQP && obj.Sample > 0 && ~isempty(obj.Memory)
                % MEX may write directly into legacy work arrays. Copy the
                % selected data before it leaves the owning controller.
                stage = {'Q','S','R','A','B','Cx','Cgx','Cgu','CgN','gx','gu', ...
                    'a','ds0','lc','uc','lb_dx','ub_dx','lb_du','ub_du'};
                condensed = {'Hc','gc','Ccx','Ccg','lcc','ucc','lxc','uxc', ...
                    'Hc_r','gc_r','Ccx_r','Ccg_r','R_block','gu_block'};
                state.qp = copyNumericFields(obj.Memory, [stage, condensed]);
                if isfield(obj.Memory, 'mem2')
                    state.qp.partial = copyNumericFields(obj.Memory.mem2, stage);
                end
            end
        end

        function reset(obj)
            input = copyNumericFields(obj.Initial, fieldnames(obj.Initial));
            mem = obj.initializeMemory();
            obj.releaseMemory();
            obj.Input = input;
            obj.Memory = mem;
            obj.Sample = 0;
            obj.Info = struct();
        end

        function delete(obj)
            for k = 1:numel(obj)
                if ~isvalid(obj(k))
                    continue;
                end
                lease = obj(k).Lease;
                obj(k).Lease = [];
                obj(k).releaseMemory();
                clear lease;
            end
        end
    end

    methods (Static)
        function assertBuildAllowed()
            if controllerCount(0) > 0
                error('MATMPC:ActiveControllers', ...
                    'Delete the live MATMPCController instances before rebuilding or clearing the MATMPC MEX files.');
            end
        end
    end

    methods (Access = private)
        function mem = initializeMemory(obj)
            mem = InitMemory(obj.Settings, obj.Options, obj.Initial);
            if any(strcmp(obj.Options.qpsolver, {'qpoases', 'qpoases_mb'}))
                mem.own_qpoases_sequence = true;
            end
        end

        function releaseMemory(obj)
            mem = obj.Memory;
            obj.Memory = [];
            if isempty(mem)
                return;
            end
            if isfield(mem, 'qp_obj')
                delete(mem.qp_obj);
            elseif isfield(mem, 'mem2') && isfield(mem.mem2, 'qp_obj')
                delete(mem.mem2.qp_obj);
            end
            % Let the scoped qpOASES cleanup close this instance's sequence.
        end
    end
end

function count = controllerCount(delta)
    persistent active
    if isempty(active)
        active = 0;
    end
    active = active + delta;
    count = active;
end

function opt = controllerOptions(passed)
    validateattributes(passed, {'struct'}, {'scalar'}, mfilename, 'options');
    opt = struct('hessian', 'Gauss_Newton', 'integrator', 'ERK4', 'qpsolver', 'qpoases', ...
        'condensing', '', 'hotstart', 'no', 'nonuniform_grid', false, 'RTI', 'yes');
    keys = fieldnames(opt);
    for i = 1:numel(keys)
        key = keys{i};
        if isfield(passed, key)
            opt.(key) = passed.(key);
        end
    end
    validateattributes(opt.nonuniform_grid, {'logical', 'double'}, {'scalar', 'binary'}, mfilename, 'nonuniform_grid');
    opt.nonuniform_grid = logical(opt.nonuniform_grid);
    choices = {{'Gauss_Newton', 'Generalized_Gauss_Newton'}, {'ERK4', 'IRK3', 'IRK3-DAE'}};
    opt.hessian = validatestring(opt.hessian, choices{1}, mfilename, 'hessian');
    opt.integrator = validatestring(opt.integrator, choices{2}, mfilename, 'integrator');
    opt.RTI = validatestring(opt.RTI, {'yes', 'no'}, mfilename, 'RTI');
    opt.hotstart = validatestring(opt.hotstart, {'yes', 'no'}, mfilename, 'hotstart');
    dense = {'qpoases', 'qpoases_mb', 'quadprog_dense', 'ipopt_dense'};
    sparse = {'osqp_sparse', 'ipopt_sparse', 'hpipm_sparse', 'hpipm_pcond'};
    partial = {'osqp_partial_sparse', 'ipopt_partial_sparse'};
    opt.qpsolver = validatestring(opt.qpsolver, [dense, sparse, partial], mfilename, 'qpsolver');
    if any(strcmp(opt.qpsolver, dense))
        condensing = {'default_full', 'blasfeo_full'};
        if strcmp(opt.qpsolver, 'qpoases_mb')
            condensing = condensing(1);
        end
    elseif any(strcmp(opt.qpsolver, sparse))
        condensing = {'no'};
    else
        condensing = {'partial_condensing'};
    end
    if isempty(opt.condensing)
        opt.condensing = condensing{1};
    else
        opt.condensing = validatestring(opt.condensing, condensing, mfilename, 'condensing');
    end
end

function settings = controllerSettings(project, requested, opt)
    validateattributes(requested, {'struct'}, {'scalar'}, mfilename, 'settings');
    addpath(fullfile(project, 'nmpc'));
    builtCore = fullfile(project, 'mex_core');
    if exist(fullfile(builtCore, ['qp_generation.', mexext]), 'file') ~= 0
        addpath(builtCore, fullfile(project, 'model_src'));
    end
    gateway = which('qp_generation');
    if isempty(gateway) || exist(gateway, 'file') ~= 3
        error('MATMPC:MissingControllerBuild', 'Build the selected model with Build_Model before creating a controller.');
    end
    buildRoot = fileparts(fileparts(gateway));
    metadata = fullfile(buildRoot, 'data', 'settings.mat');
    if exist(metadata, 'file') ~= 2
        error('MATMPC:MissingControllerBuild', 'The active MATMPC MEX directory needs its generated data/settings.mat.');
    end
    model = load(metadata, 'settings');
    geometry = {'model','nx','nu','nz','np','ny','nyN','nc','ncN','nbx','nbu','nbx_idx','nbu_idx'};
    for i = 1:numel(geometry)
        key = geometry{i};
        if isfield(requested, key) && ~isequal(requested.(key), model.settings.(key))
            error('MATMPC:ModelSettingsMismatch', 'The active model build has a different settings.%s.', key);
        end
    end
    settings = copyFields(struct(), model.settings, [geometry, {'Ts_st'}]);
    settings.N = 80;
    settings.r = 10;
    settings = copyFields(settings, requested, {'Ts_st', 'N', 'r', 'index_T'});
    validateattributes(settings.Ts_st, {'double'}, {'real', 'finite', 'scalar', 'positive'}, mfilename, 'Ts_st');
    validateattributes(settings.N, {'double'}, {'real', 'finite', 'scalar', 'integer', 'positive'}, mfilename, 'N');
    validateattributes(settings.r, {'double'}, {'real', 'finite', 'scalar', 'integer', 'positive'}, mfilename, 'r');
    if mod(settings.N, 5) == 0
        settings.N2 = settings.N / 5;
    else
        settings.N2 = settings.N;
    end
    settings = copyFields(settings, requested, {'N2'});
    validateattributes(settings.N2, {'double'}, {'real', 'finite', 'scalar', 'integer', 'positive'}, mfilename, 'N2');
    if isfield(settings, 'index_T')
        validateattributes(settings.index_T, {'double'}, {'real', 'finite', 'vector', 'integer'}, mfilename, 'index_T');
    end
    addSolverPath(buildRoot, opt.qpsolver);
end

function addSolverPath(root, name)
    if ispc
        platform = 'win64';
    elseif ismac
        platform = 'mac';
    else
        platform = 'linux';
    end
    if any(strcmp(name, {'qpoases', 'qpoases_mb'}))
        solver = 'qpoases';
        mexName = 'qpOASES_sequence';
    elseif any(strcmp(name, {'osqp_sparse', 'osqp_partial_sparse'}))
        solver = 'osqp';
        mexName = 'osqp_mex';
    else
        return;
    end
    folder = fullfile(root, 'solver', platform, solver);
    if exist(fullfile(folder, [mexName, '.', mexext]), 'file') ~= 0
        addpath(folder, '-end');
    end
    if exist(mexName, 'file') ~= 3
        error('MATMPC:MissingControllerSolver', 'Add %s built for %s to the MATLAB path.', mexName, mexext);
    end
end

function input = validateInput(s, data)
    validateattributes(data, {'struct'}, {'scalar'}, mfilename, 'input');
    names = [{'x0'}, iterateNames(), tuningNames()];
    for i = 1:numel(names)
        name = names{i};
        if ~isfield(data, name)
            error('MATMPC:MissingControllerInput', 'Field input.%s is required.', name);
        end
        input.(name) = numericCopy(checkInputField(s, name, data.(name)));
    end
    checkBounds(input);
end

function data = applyUpdates(s, data, updates)
    validateattributes(updates, {'struct'}, {'scalar'}, mfilename, 'updates');
    keys = fieldnames(updates);
    allowed = tuningNames();
    for i = 1:numel(keys)
        key = keys{i};
        if ~ismember(key, allowed)
            error('MATMPC:UnsupportedControllerUpdate', 'Field updates.%s is not runtime tuning data.', key);
        end
        data.(key) = checkInputField(s, key, updates.(key));
    end
    if ~isempty(keys)
        checkBounds(data);
    end
end

function value = checkInputField(s, name, value)
    ncon = s.N * s.nc + s.ncN;
    shapes = struct('x0', [s.nx, 1], 'x', [s.nx, s.N+1], 'u', [s.nu, s.N], ...
        'z', [s.nz, s.N], 'lambda', [s.nx, s.N+1], 'mu', [ncon, 1], ...
        'mu_x', [s.N * s.nbx, 1], 'mu_u', [s.N * s.nu, 1], ...
        'y', [s.ny, s.N], 'yN', [s.nyN, 1], 'od', [max(s.np, 1), s.N+1], ...
        'W', [s.ny, s.N], 'WN', [s.nyN, 1], 'lbu', [s.nu, s.N], 'ubu', [s.nu, s.N], ...
        'lbx', [s.nbx, s.N], 'ubx', [s.nbx, s.N], 'lb', [ncon, 1], 'ub', [ncon, 1]);
    dimensions = shapes.(name);
    attributes = {'real', 'nonsparse', 'nonnan'};
    if ~any(strcmp(name, {'lbu', 'ubu', 'lbx', 'ubx', 'lb', 'ub'}))
        attributes{end+1} = 'finite';
    end
    if prod(dimensions) == 0 && isempty(value)
        value = reshape(value, dimensions);
    end
    validateattributes(value, {'double'}, [attributes, {'size', dimensions}], mfilename, name);
end

function checkBounds(input)
    pairs = {'lbu','ubu'; 'lbx','ubx'; 'lb','ub'};
    for i = 1:size(pairs, 1)
        lower = input.(pairs{i, 1});
        upper = input.(pairs{i, 2});
        if any(lower(:) > upper(:) | lower(:) == Inf | upper(:) == -Inf)
            error('MATMPC:InvalidControllerBounds', '%s and %s do not describe valid bounds.', pairs{i, :});
        end
    end
end

function names = iterateNames()
    names = {'x','u','z','lambda','mu','mu_x','mu_u'};
end

function names = tuningNames()
    names = {'y','yN','od','W','WN','lbu','ubu','lbx','ubx','lb','ub'};
end

function result = copyFields(result, source, names)
    for i = 1:numel(names)
        name = names{i};
        if isfield(source, name)
            result.(name) = source.(name);
        end
    end
end

function result = copyNumericFields(source, names)
    result = struct();
    for i = 1:numel(names)
        name = names{i};
        if isfield(source, name)
            result.(name) = numericCopy(source.(name));
        end
    end
end

function result = numericCopy(value)
    result = zeros(size(value), 'like', value);
    % Separate assignments materialize the result's MATLAB-owned storage.
    if ~isempty(value)
        result(1) = value(1);
        result(2:end) = value(2:end);
    end
end

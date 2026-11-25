function smoke_inverted_pendulum()
%SMOKE_INVERTED_PENDULUM Build and run the default model for five samples.

projectRoot = fileparts(fileparts(mfilename('fullpath')));
originalPath = path;
restorePath = onCleanup(@() path(originalPath));
addpath(projectRoot);

solverRoot = fullfile(projectRoot, 'solver');
if ismac
    solverDir = 'mac';
elseif isunix
    solverDir = 'linux';
elseif ispc
    solverDir = 'win64';
else
    error('MATMPC:UnsupportedPlatform', 'No bundled qpOASES build for this platform.');
end
addpath(fullfile(solverRoot, solverDir, 'qpoases'), '-end');
assert(exist('qpOASES_sequence', 'file') == 3, ...
    'A qpOASES_sequence MEX file matching %s is required.', mexext);

settings = Build_Model('InvertedPendulum');
addpath(fullfile(projectRoot, 'mex_core'), fullfile(projectRoot, 'nmpc'));

settings.N = 80;
settings.N2 = settings.N / 5;
settings.r = 10;
opt.hessian = 'Gauss_Newton';
opt.integrator = 'ERK4';
opt.condensing = 'default_full';
opt.qpsolver = 'qpoases';
opt.hotstart = 'no';
opt.nonuniform_grid = 0;
opt.RTI = 'yes';

[input, data] = InitData(settings);
input.y = repmat(data.REF', 1, settings.N);
input.yN = data.REF(1:settings.nyN)';
mem = InitMemory(settings, opt, input);
% The MEX solver retains the active-set sequence between successive calls.
cleanupSolver = onCleanup(@() clearQpOases());
state = input.x0;
boundTolerance = 1e-6;

for step = 1:5
    mem.iter = step;
    input.x0 = state;
    [output, mem] = mpc_nmpcsolver(input, settings, mem, opt);
    control = output.u(:, 1);
    [state, ~] = Simulate_System(state, control, output.z(:, 1), input.od(:, 1), mem, settings);
    assert(all(isfinite([state(:); control(:); output.info.OptCrit])), ...
        'Nonfinite state, control, or solver residual at step %d.', step);
    assert(output.info.iteration_num == 1, 'Expected one RTI iteration at step %d.', step);
    assert(all(control >= input.lbu(:, 1) - boundTolerance & ...
        control <= input.ubu(:, 1) + boundTolerance), 'Control bound violated at step %d.', step);
    boundedState = state(settings.nbx_idx);
    assert(all(boundedState >= input.lbx(:, 1) - boundTolerance & ...
        boundedState <= input.ubx(:, 1) + boundTolerance), 'State bound violated at step %d.', step);

    input.x = output.x;
    input.u = output.u;
    input.z = output.z;
    input.lambda = output.lambda;
    input.mu = output.mu;
    input.mu_x = output.mu_x;
    input.mu_u = output.mu_u;
end

qpOASES_sequence('c', mem.warm_start);
fprintf('InvertedPendulum smoke test passed 5 closed-loop steps.\n');
end

function clearQpOases()
clear qpOASES_sequence;
end

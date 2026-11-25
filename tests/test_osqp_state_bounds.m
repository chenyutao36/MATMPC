function tests = test_osqp_state_bounds
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    testCase.addTeardown(@path, path);
    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(root, fullfile(root, 'nmpc'));
    if ispc
        platform = 'win64';
    elseif ismac
        platform = 'mac';
    else
        platform = 'linux';
    end
    solverDir = fullfile(root,'solver',platform,'osqp');
    available = exist(fullfile(solverDir,['osqp_mex.' mexext]),'file') ~= 0;
    if available
        addpath(solverDir);
    end
    testCase.TestData.hasOSQP = available;
end

function testSparseStateBounds(testCase)
    assumeTrue(testCase,testCase.TestData.hasOSQP,'Bundled OSQP MEX is unavailable on this platform.');
    settings = makeSettings();
    mem = InitMemory(settings,makeOptions('osqp_sparse','no'),struct());
    verifyFalse(testCase, isfield(mem, 'Hc'));
    checkIndependentBounds(testCase, mem, settings);
end

function testPartialSparseStateBounds(testCase)
    assumeTrue(testCase,testCase.TestData.hasOSQP,'Bundled OSQP MEX is unavailable on this platform.');
    settings = makeSettings();
    mem = InitMemory(settings,makeOptions('osqp_partial_sparse','partial_condensing'),struct());
    verifyFalse(testCase, isfield(mem, 'Hc'));
    checkIndependentBounds(testCase, mem.mem2, mem.settings2);
end

function testInvalidPartialHorizon(testCase)
    settings = makeSettings();
    opt = makeOptions('osqp_partial_sparse','partial_condensing');
    for invalid = [0, 1.5, 3]
        settings.N2 = invalid;
        verifyError(testCase,@() InitMemory(settings,opt,struct()),'MATMPC:InvalidCondensingHorizon');
    end
end

function settings = makeSettings()
    settings = struct('Ts_st',0.05,'nx',4,'nu',1,'nz',0,'ny',4,'nyN',4,'np',0, ...
        'nbx',2,'nbu',1,'nbx_idx',[2,4],'nbu_idx',1,'nc',0,'ncN',0,'N',4,'N2',2,'r',1);
end

function opt = makeOptions(solver,condensing)
    opt = struct('hotstart','no','condensing',condensing,'qpsolver',solver, ...
        'hessian','Gauss_Newton','integrator','ERK4','RTI','yes','nonuniform_grid',0);
end

function checkIndependentBounds(testCase, mem, settings)
    [~, ~, pattern] = sparse_qp_layout(settings, true);
    horizon = settings.N;
    stateRows = (horizon + 1) * settings.nx + horizon * settings.nu + (1:horizon * settings.nbx);
    A = pattern(stateRows, :);
    cols = bsxfun(@plus, settings.nbx_idx(:), settings.nx * (1:horizon));
    expectedA = sparse((1:numel(cols))',cols(:),ones(numel(cols),1),size(A,1),size(A,2));
    verifyEqual(testCase,A,expectedA);

    % The bounds of the two states are disjoint at every stage.
    % The optimizer of 0.5*x'*x is the projection of zero onto each interval.
    lower = repmat([-2;0.25],horizon,1);
    upper = repmat([-1;0.5],horizon,1);
    solver = osqp;
    solver.setup(speye(size(A,2)),zeros(size(A,2),1),A,lower,upper,mem.osqp_options);
    sol = solver.solve();
    assertEqual(testCase,sol.info.status,'solved');
    verifyEqual(testCase,A*sol.x,repmat([-1;0.25],horizon,1),'AbsTol',1e-4);
end

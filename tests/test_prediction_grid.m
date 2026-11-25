function tests = test_prediction_grid
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    testCase.addTeardown(@path, path);
    root = fileparts(fileparts(mfilename('fullpath')));
    if ismac
        platform = 'mac';
    elseif ispc
        platform = 'win64';
    else
        platform = 'linux';
    end
    addpath(root, fullfile(root,'solver',platform,'qpoases'));
end

function testNonuniformGridUsesPredictionIntervals(testCase)
    settings = exampleSettings();
    settings.N = 3;
    settings.r = 3;
    settings.index_T = [0; 1; 4; 6];
    mem = InitMemory(settings,exampleOptions(),struct());
    verifyTrue(testCase,mem.nonuniform_grid);
    verifyEqual(testCase,mem.index_T,[0,1,4,6]);
    verifyFalse(testCase,isfield(mem,'T'));
end

function testNonuniformGridRejectsInvalidBoundaries(testCase)
    settings = exampleSettings();
    settings.N = 3;
    settings.r = 3;
    for boundaries = {[0 1 2], [1 2 3 4], [0 1 1 3], [0 1 1.5 3], [0 1 2 Inf]}
        settings.index_T = boundaries{1};
        verifyError(testCase,@() InitMemory(settings,exampleOptions(),struct()),'MATMPC:InvalidGrid');
    end
    settings.index_T = [0 1 2 3];
    settings.N = 4;
    verifyError(testCase,@() InitMemory(settings,exampleOptions(),struct()),'MATMPC:InvalidGrid');
end

function testMoveBlockingUsesSeparateControlBlocks(testCase)
    settings = exampleSettings();
    settings.index_T = [0 2 5];
    opt = exampleOptions();
    opt.nonuniform_grid = false;
    opt.qpsolver = 'qpoases_mb';
    mem = InitMemory(settings,opt,struct());
    verifyFalse(testCase,mem.nonuniform_grid);
    verifyEqual(testCase,mem.T,kron([1 0; 1 0; 0 1; 0 1; 0 1],eye(settings.nu)));
    verifySize(testCase,mem.Hc_r,[settings.r*settings.nu,settings.r*settings.nu]);
    verifyFalse(testCase,any(isfield(mem,{'Hc','gc','Ccx','Ccg','Gr','L'})));
    settings.index_T = [0 2 4];
    verifyError(testCase,@() InitMemory(settings,opt,struct()),'MATMPC:InvalidGrid');
    opt.nonuniform_grid = true;
    verifyError(testCase,@() InitMemory(settings,opt,struct()),'MATMPC:IncompatibleGrid');
end

function testLevenbergMarquardtOption(testCase)
    settings = exampleSettings();
    opt = exampleOptions();
    opt.nonuniform_grid = false;
    mem = InitMemory(settings,opt,struct());
    verifyEqual(testCase,mem.levenberg_marquardt,0);
    opt.levenberg_marquardt = single(0.25);
    mem = InitMemory(settings,opt,struct());
    verifyEqual(testCase,mem.levenberg_marquardt,0.25);
    for invalid = {-1, NaN, Inf, [0 1], 1i, sparse(1), '1', []}
        opt.levenberg_marquardt = invalid{1};
        verifyError(testCase,@() InitMemory(settings,opt,struct()),'MATMPC:InvalidLevenbergMarquardt');
    end
end

function settings = exampleSettings()
    settings = struct('Ts_st',0.05,'nx',1,'nu',2,'nz',0,'ny',1,'nyN',1,'np',0, ...
        'nbx',0,'nbu',0,'nbx_idx',[],'nbu_idx',[],'nc',0,'ncN',0,'N',5,'N2',1,'r',2);
end

function opt = exampleOptions()
    opt = struct('hotstart','no','condensing','default_full','qpsolver','quadprog_dense', ...
        'hessian','Gauss_Newton','integrator','ERK4','RTI','yes','nonuniform_grid',true);
end

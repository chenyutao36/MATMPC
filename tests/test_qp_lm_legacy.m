function tests = test_qp_lm_legacy
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.assumeNotEmpty(mex.getCompilerConfigurations('C', 'Selected'), ...
    'Configure a MATLAB C compiler with mex -setup C to run the legacy LM tests.');
testDir = fileparts(mfilename('fullpath'));
root = fileparts(testDir);
mexCore = fullfile(root, 'mex_core');
portableCore = fullfile(root, 'core');
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root));
if ismac
    platform = 'mac';
elseif ispc
    platform = 'win64';
else
    platform = 'linux';
end
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'solver', platform, 'qpoases')));
tmp = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(tmp.Folder));
testCase.addTeardown(@unloadFixtures);
mexSupport = fullfile(mexCore, {'sim.c', 'erk.c', 'irk_mex.c', 'irk_ode.c', 'irk_dae.c', 'mpc_common.c'});
coreSupport = fullfile(portableCore, {'qp_gauss_newton.c', 'erk_kernel.c', 'irk_kernel.c', 'matmpc_lu.c', 'matmpc_la.c'});
for target = {'qp_generation', 'qp_generation_ngrid', 'qp_generation_mb'}
    name = target{1};
    mex('-outdir', tmp.Folder, ['-I' mexCore], ['-I' portableCore], ...
        fullfile(mexCore, [name '.c']), mexSupport{:}, coreSupport{:}, ...
        fullfile(testDir, 'fixtures', 'scalar_ode.c'), '-lmwblas', '-lmwlapack');
end
end

function unloadFixtures
clear qp_generation qp_generation_ngrid qp_generation_mb
end

function testAdditiveShiftAfterFloorAcrossLegacyModes(testCase)
shift = 0.25;
for entry = {'qp_generation', 'qp_generation_ngrid', 'qp_generation_mb'}
    for hessian = {'Gauss_Newton', 'Generalized_Gauss_Newton'}
        for integrator = {'ERK4', 'IRK3', 'IRK3-DAE'}
            [input, settings, opt] = problem(entry{1}, hessian{1}, integrator{1});
            legacy = InitMemory(settings, opt, input);
            withoutShift = InitMemory(settings, opt, input);
            legacy.reg = 4;
            withoutShift.reg = 4;
            if isfield(withoutShift, 'levenberg_marquardt')
                withoutShift = rmfield(withoutShift, 'levenberg_marquardt');
            end
            legacy.levenberg_marquardt = shift;
            feval(entry{1}, input, settings, withoutShift);
            feval(entry{1}, input, settings, legacy);

            testCase.verifyEqual(legacy.Q, [4, 4, 4, 8] + shift, 'AbsTol', 1e-12);
            if strcmp(entry{1}, 'qp_generation_mb')
                testCase.verifyEqual(withoutShift.R, [5, 2, 3], 'AbsTol', 1e-12);
                testCase.verifyEqual(legacy.R, [5, 2, 3] + shift, 'AbsTol', 1e-12);
                testCase.verifyEqual(legacy.R_block, [7, 4] + [2, 1] * shift, 'AbsTol', 1e-12);
                testCase.verifyEqual(legacy.gu_block, withoutShift.gu_block);
            else
                testCase.verifyEqual(withoutShift.R, [5, 4, 4], 'AbsTol', 1e-12);
                testCase.verifyEqual(legacy.R, [5, 4, 4] + shift, 'AbsTol', 1e-12);
            end
            testCase.verifyEqual(legacy.S, withoutShift.S);
            testCase.verifyEqual(legacy.gx, withoutShift.gx);
            testCase.verifyEqual(legacy.gu, withoutShift.gu);
            testCase.verifyEqual(legacy.A, withoutShift.A);
            testCase.verifyEqual(legacy.B, withoutShift.B);
        end
    end
end
end

function testInvalidShiftDoesNotOverwriteWorkspace(testCase)
identifier = 'MATMPC:QPGeneration:LevenbergMarquardt';
for entry = {'qp_generation', 'qp_generation_ngrid', 'qp_generation_mb'}
    [input, settings, opt] = problem(entry{1}, 'Gauss_Newton', 'IRK3');
    mem = InitMemory(settings, opt, input);
    mem.A = mem.A + 42;
    for value = {-1, NaN, Inf, [0, 1], sparse(1), 1 + 1i, '0'}
        mem.levenberg_marquardt = value{1};
        testCase.verifyError(@() runLegacy(entry{1}, input, settings, mem), identifier);
        testCase.verifyEqual(mem.A, 42 * ones(size(mem.A)));
    end
    if strcmp(entry{1}, 'qp_generation_mb')
        mem.levenberg_marquardt = realmax;
        testCase.verifyError(@() runLegacy(entry{1}, input, settings, mem), ...
            'MATMPC:QPGeneration:Configuration');
        testCase.verifyEqual(mem.A, 42 * ones(size(mem.A)));
    end
    mem.levenberg_marquardt = 0;
    feval(entry{1}, input, settings, mem);
    testCase.verifyEqual(mem.A, ones(size(mem.A)), 'AbsTol', 1e-12);
end
end

function testShiftRejectsNonfiniteResult(testCase)
[input, settings, opt] = problem('qp_generation', 'Gauss_Newton', 'ERK4');
input.W(1) = realmax;
mem = InitMemory(settings, opt, input);
mem.levenberg_marquardt = realmax;
testCase.verifyError(@() runLegacy('qp_generation', input, settings, mem), ...
    'MATMPC:QPGeneration:Assembly');
end

function runLegacy(entry, input, settings, mem)
feval(entry, input, settings, mem);
end

function [input, settings, opt] = problem(entry, hessian, integrator)
settings = struct('Ts_st', 0.25, 'nx', 1, 'nu', 1, 'nz', double(strcmp(integrator, 'IRK3-DAE')), ...
    'ny', 1, 'nyN', 1, 'np', 0, 'nbx', 0, 'nbu', 0, 'nbx_idx', [], 'nbu_idx', [], ...
    'nc', 0, 'ncN', 0, 'N', 3, 'N2', 1, 'r', 3, 'index_T', [0, 1, 4, 9]);
opt = struct('hotstart', 'no', 'condensing', 'default_full', 'qpsolver', 'quadprog_dense', ...
    'hessian', hessian, 'integrator', integrator, 'RTI', 'yes', ...
    'nonuniform_grid', strcmp(entry, 'qp_generation_ngrid'));
if strcmp(entry, 'qp_generation_mb')
    opt.qpsolver = 'qpoases_mb';
    settings.r = 2;
    settings.index_T = [0, 2, 3];
end
input = struct('x0', 1, 'x', [1, 2, 4, 5], 'u', [1, 2, 1], ...
    'z', zeros(settings.nz, settings.N), 'od', zeros(1, settings.N+1), ...
    'y', zeros(1, settings.N), 'yN', 1, 'W', [5, 2, 3], 'WN', 8, ...
    'lbu', -Inf(1, settings.N), 'ubu', Inf(1, settings.N), ...
    'lb', [], 'ub', [], 'lbx', [], 'ubx', []);
end

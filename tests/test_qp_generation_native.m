function tests = test_qp_generation_native
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.assumeNotEmpty(mex.getCompilerConfigurations('C', 'Selected'), ...
    'Configure a MATLAB C compiler with mex -setup C to run QP assembly tests.');
testDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(testDir);
core = fullfile(projectRoot, 'core');
mexCore = fullfile(projectRoot, 'mex_core');
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(projectRoot));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(projectRoot, 'nmpc')));
tmp = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(tmp.Folder));
testCase.addTeardown(@unloadFixture);
sources = fullfile(mexCore, {'sim.c', 'erk.c', 'irk_ode.c', 'irk_dae.c', 'mpc_common.c'});
portable = fullfile(core, {'qp_gauss_newton.c', 'erk_kernel.c', 'matmpc_la.c'});
for entry = {'qp_generation', 'qp_generation_ngrid'}
    name = entry{1};
    mex('-outdir', tmp.Folder, ['-I' core], ['-I' mexCore], '-output', name, ...
        fullfile(mexCore, [name '.c']), sources{:}, portable{:}, ...
        fullfile(testDir, 'fixtures', 'scalar_ode.c'), '-lmwblas', '-lmwlapack');
end
end

function unloadFixture
clear qp_generation qp_generation_ngrid
end

function testExplicitOutputsMatchLegacyBlockByBlock(testCase)
[input, settings, opt] = scalarProblem();
mem = InitMemory(settings, opt, input);
legacy = InitMemory(settings, opt, input);
qp = qp_generation(input, settings, mem);
qp_generation(input, settings, legacy);
fields = fieldnames(qp);
testCase.verifyNumElements(fields, 18);
for i = 1:numel(fields)
    name = fields{i};
    testCase.verifyEqual(qp.(name), legacy.(name), 'AbsTol', 1e-12, name);
end
testCase.verifyEqual(qp.A, ones(1, settings.N), 'AbsTol', 1e-12);
testCase.verifyEqual(qp.B, settings.Ts_st * ones(1, settings.N), 'AbsTol', 1e-12);
testCase.verifyEqual(qp.lb_dx, input.lbx(:) - input.x(2:end).', 'AbsTol', 1e-12);
testCase.verifyEqual(mem.A, zeros(size(mem.A)));
testCase.verifyEqual(mem.R, zeros(size(mem.R)));
testCase.verifyEqual(input.x, [1, 2, 4, 5]);

updated = prepare_uniform_qp(input, settings, mem);
for i = 1:numel(fields)
    name = fields{i};
    testCase.verifyEqual(updated.(name), qp.(name), 'AbsTol', 1e-12, name);
end
testCase.verifyEqual(mem.A, zeros(size(mem.A)));
testCase.verifyEqual(updated.warm_start, mem.warm_start);
end

function testNonuniformOutputsMatchLegacyAndKeepThePredictionHorizon(testCase)
[input, settings, opt] = scalarProblem();
settings.index_T = [0, 1, 4, 9];
opt.nonuniform_grid = true;
mem = InitMemory(settings, opt, input);
legacy = InitMemory(settings, opt, input);
qp = qp_generation(input, settings, mem);
qp_generation_ngrid(input, settings, legacy);
updated = prepare_shooting_qp(input, settings, mem);
fields = fieldnames(qp);
for i = 1:numel(fields)
    name = fields{i};
    testCase.verifyEqual(qp.(name), legacy.(name), 'AbsTol', 1e-12, name);
    testCase.verifyEqual(updated.(name), qp.(name), 'AbsTol', 1e-12, name);
end

duration = settings.Ts_st * diff(settings.index_T);
testCase.verifyEqual(qp.B, duration, 'AbsTol', 1e-12);
testCase.verifyEqual(qp.a, input.x(1:end-1) + input.u .* duration - input.x(2:end), 'AbsTol', 1e-12);
testCase.verifyEqual(qp.lb_dx, input.lbx(:) - input.x(2:end).', 'AbsTol', 1e-12);
testCase.verifyEqual(qp.ub_dx, input.ubx(:) - input.x(2:end).', 'AbsTol', 1e-12);
testCase.verifySize(qp.A, [1, settings.N]);
testCase.verifySize(qp.gx, [1, settings.N + 1]);
testCase.verifyEqual(mem.A, zeros(size(mem.A)));

invalid = mem;
invalid.index_T = [0, 1, 1, 9];
testCase.verifyError(@() explicitQp(input, settings, invalid), 'MATMPC:Integrator:Grid');
mem.nonuniform_grid = false;
uniform = prepare_shooting_qp(input, settings, mem);
testCase.verifyEqual(uniform.B, settings.Ts_st * ones(1, settings.N), 'AbsTol', 1e-12);
end

function testInvalidInputsLeaveDefaultQpUsable(testCase)
[input, settings, opt] = scalarProblem();
mem = InitMemory(settings, opt, input);
bad = input;
bad.u = single(bad.u);
testCase.verifyError(@() explicitQp(bad, settings, mem), 'MATMPC:QPGeneration:Input');
bad = input;
bad.x = bad.x(:, 2:end);
testCase.verifyError(@() explicitQp(bad, settings, mem), 'MATMPC:QPGeneration:Input');
badSettings = settings;
badSettings.nbx_idx = 2;
testCase.verifyError(@() explicitQp(input, badSettings, mem), 'MATMPC:QPGeneration:StateBounds');
qp = qp_generation(input, settings, mem);
testCase.verifyEqual(qp.B, settings.Ts_st * ones(1, settings.N), 'AbsTol', 1e-12);
end

function testGgnAndIrkKeepTheirExistingPaths(testCase)
[input, settings, opt] = scalarProblem();
for selected = {'Generalized_Gauss_Newton', 'IRK3'}
    options = opt;
    if strcmp(selected{1}, 'IRK3')
        options.integrator = 'IRK3';
    else
        options.hessian = selected{1};
    end
    mem = InitMemory(settings, options, input);
    testCase.verifyError(@() explicitQp(input, settings, mem), 'MATMPC:QPGeneration:Mode');
    updated = prepare_uniform_qp(input, settings, mem);
    testCase.verifyEqual(updated.A, ones(1, settings.N), 'AbsTol', 1e-12);
    testCase.verifyEqual(updated.B, settings.Ts_st * ones(1, settings.N), 'AbsTol', 1e-12);
    testCase.verifyEqual(updated.R, input.W, 'AbsTol', 1e-12);
end
end

function testNonuniformGgnAndIrkKeepTheirExistingPaths(testCase)
[input, settings, opt] = scalarProblem();
settings.index_T = [0, 1, 4, 9];
opt.nonuniform_grid = true;
duration = settings.Ts_st * diff(settings.index_T);
for selected = {'Generalized_Gauss_Newton', 'IRK3'}
    options = opt;
    if strcmp(selected{1}, 'IRK3')
        options.integrator = 'IRK3';
    else
        options.hessian = selected{1};
    end
    mem = InitMemory(settings, options, input);
    testCase.verifyError(@() explicitQp(input, settings, mem), 'MATMPC:QPGeneration:Mode');
    updated = prepare_shooting_qp(input, settings, mem);
    testCase.verifyEqual(updated.A, ones(1, settings.N), 'AbsTol', 1e-12);
    testCase.verifyEqual(updated.B, duration, 'AbsTol', 1e-12);
    testCase.verifyEqual(updated.R, input.W, 'AbsTol', 1e-12);
end
end

function explicitQp(input, settings, mem)
unused = qp_generation(input, settings, mem); %#ok<NASGU>
end

function [input, settings, opt] = scalarProblem()
settings = struct('Ts_st', 0.25, 'nx', 1, 'nu', 1, 'nz', 0, 'np', 0, 'ny', 1, 'nyN', 1, ...
    'nbx', 1, 'nbu', 0, 'nbx_idx', 1, 'nbu_idx', [], 'nc', 0, 'ncN', 0, 'N', 3, 'N2', 1, 'r', 3);
opt = struct('hotstart', 'no', 'condensing', 'default_full', 'qpsolver', 'quadprog_dense', ...
    'hessian', 'Gauss_Newton', 'integrator', 'ERK4', 'RTI', 'yes', 'nonuniform_grid', false);
input = struct('x0', 0.5, 'x', [1, 2, 4, 5], 'u', [2, -1, 0.5], 'od', zeros(1, 4), 'z', zeros(0, 3), ...
    'y', [0, 2, 3], 'yN', 3, 'W', [5, 2, 3], 'WN', 4, 'lbu', [-2, -3, -4], 'ubu', [4, 5, 6], ...
    'lbx', [-1, -2, -3], 'ubx', [2, 3, 4], 'lb', [], 'ub', []);
end

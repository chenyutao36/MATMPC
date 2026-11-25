function tests = test_sqp_stage_consistency
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.assumeNotEmpty(mex.getCompilerConfigurations('C', 'Selected'), ...
    'Configure a MATLAB C compiler with mex -setup C to run the SQP MEX tests.');

testDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(testDir);
coreDir = fullfile(projectRoot, 'mex_core');
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(projectRoot));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(projectRoot, 'nmpc')));
tmp = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(tmp.Folder));
testCase.addTeardown(@unloadFixtures);
options = {'-outdir', tmp.Folder, ['-I' coreDir]};
support = fullfile(coreDir, {'sim.c', 'erk.c', 'irk_ode.c', 'irk_dae.c', 'mpc_common.c'});
model = fullfile(testDir, 'fixtures', 'scalar_ode.c');
for target = {'qp_generation_ngrid', 'Line_search', 'solution_info'}
    gateway = target{1};
    mex(options{:}, '-output', gateway, fullfile(coreDir, [gateway '.c']), ...
        support{:}, model, '-lmwblas', '-lmwlapack');
end
end

function unloadFixtures
clear qp_generation_ngrid Line_search solution_info
end

function testUnequalGridMatchesAnalyticDynamicsAndSensitivity(testCase)
for integrator = {'ERK4', 'IRK3', 'IRK3-DAE'}
    [input, settings, mem, durations] = scalarProblem(integrator{1});
    input.y = input.u + durations;
    input.WN = 1;
    input.yN = input.x(end) - 1;
    input.lambda = [-1, 1, 1, 1];
    qp_generation_ngrid(input, settings, mem);
    testCase.verifyEqual(mem.a, zeros(1, settings.N), 'AbsTol', 1e-12);
    testCase.verifyEqual(mem.A, ones(1, settings.N), 'AbsTol', 1e-12);
    testCase.verifyEqual(mem.B, durations, 'AbsTol', 1e-12);

    [eq_res, ineq_res, KKT, objective, z] = solution_info(input, settings, mem);
    testCase.verifyEqual([eq_res, ineq_res, KKT], [0, 0, 0], 'AbsTol', 1e-12);
    testCase.verifyTrue(isfinite(objective));
    testCase.verifyEqual(z, zeros(settings.nz, settings.N), 'AbsTol', 1e-12);
end
end

function testLineSearchAcceptsFeasibleUnequalGridStep(testCase)
[input, settings, mem, durations] = scalarProblem('ERK4');
input.x = zeros(1, settings.N + 1);
input.u = zeros(1, settings.N);
input.y = ones(1, settings.N);
mem.dx = [0, cumsum(durations)];
mem.du = ones(1, settings.N);
mem.mu_merit = 100;
mem.sqp_it = 1;
qp_generation_ngrid(input, settings, mem);
[alpha, merit] = Line_search(mem, input, settings);
testCase.verifyEqual(alpha, 1);
testCase.verifyEqual(merit, 100);
end

function testUniformDiagnosticsIgnoreMoveBlockBoundaries(testCase)
[input, settings, mem] = scalarProblem('ERK4');
mem.nonuniform_grid = false;
mem.index_T = [0, 2, settings.N];
input.x = settings.Ts_st * (0:settings.N);
input.z = zeros(1, settings.N); % Empty-state placeholder in the Simulink bus.
mem.sqp_it = 0;
[updated, newMem, eq_res] = mpc_sqp_step(input, settings, mem);
testCase.verifyEqual(eq_res, 0, 'AbsTol', 1e-12);
testCase.verifyEqual(newMem.alpha, 1);
testCase.verifySize(updated.z, size(input.z));

input.x0 = input.x(1) + 1;
[eq_res, ~, ~, ~] = solution_info(input, settings, mem);
testCase.verifyEqual(eq_res, 1, 'AbsTol', 1e-12);
end

function testAllRejectedPreservesPrimalDualAndAlgebraicState(testCase)
[input, settings, mem, durations] = scalarProblem('IRK3-DAE');
mem.num_stages = 1;
mem.A_tab = 1;
mem.B_tab = 1;
mem.JFK = [mem.h, 0];
input.z = 3 + (1:settings.N);
input.lambda = [-1, 1, 1, 1];
input.mu_u = (1:settings.N)';
input.mu_x = 0;
input.mu = 0;
mem.dx = [0; cumsum(durations(:))];
mem.du = ones(settings.N, 1);
mem.lambda_new = 5 + zeros(settings.N + 1, 1);
mem.mu_u_new = 6 + zeros(settings.N, 1);
mem.q_dual = 7 + zeros(1, settings.N + 1);
mem.dmu = 8 + zeros(settings.N, 1);
mem.mu_merit = 2;
mem.sqp_it = 1;
qp_generation_ngrid(input, settings, mem);

[rejected, rejectedMem, eq_res, ineq_res, KKT, objective] = mpc_sqp_step(input, settings, mem);
testCase.verifyEqual(rejectedMem.alpha, 0);
testCase.verifyEqual(rejected, input);
testCase.verifyEqual(rejected.x, [0, cumsum(durations)]);
testCase.verifyEqual(rejected.u, ones(1, settings.N));
testCase.verifyEqual(rejected.z, 3 + (1:settings.N));
testCase.verifyEqual(rejected.lambda, [-1, 1, 1, 1]);
testCase.verifyEqual(rejected.mu_u, (1:settings.N)');
testCase.verifyEqual(rejectedMem.q_dual, 7 + zeros(1, settings.N + 1));
testCase.verifyEqual(rejectedMem.dmu, 8 + zeros(settings.N, 1));
testCase.verifyEqual(rejectedMem.mu_merit, 2);
testCase.verifyTrue(all(isfinite([eq_res, ineq_res, KKT, objective])));

% The same uphill step remains a full step when running in RTI mode.
mem.sqp_maxit = 1;
[accepted, acceptedMem] = mpc_sqp_step(input, settings, mem);
testCase.verifyEqual(acceptedMem.alpha, 1);
testCase.verifyEqual(accepted.u, input.u + 1);
testCase.verifyEqual(accepted.z, zeros(settings.nz, settings.N), 'AbsTol', 1e-12);
testCase.verifyEqual(accepted.lambda, mem.lambda_new.');
testCase.verifyEqual(accepted.mu_u, mem.mu_u_new);
testCase.verifyEqual([accepted.mu_x, accepted.mu], [0, 0]);
testCase.verifyEqual(input.z, 3 + (1:settings.N));
end

function [input, settings, mem, durations] = scalarProblem(integrator)
settings = struct('Ts_st', 0.25, 'nx', 1, 'nu', 1, 'nz', double(strcmp(integrator, 'IRK3-DAE')), ...
    'ny', 1, 'nyN', 1, 'np', 0, 'nbx', 0, 'nbu', 0, 'nbx_idx', [], 'nbu_idx', [], ...
    'nc', 0, 'ncN', 0, 'N', 3, 'N2', 1, 'r', 3, 'index_T', [0, 1, 4, 6]);
opt = struct('hotstart', 'no', 'condensing', 'default_full', 'qpsolver', 'quadprog_dense', ...
    'hessian', 'Gauss_Newton', 'integrator', integrator, 'RTI', 'no', 'nonuniform_grid', true);
mem = InitMemory(settings, opt, struct());
durations = settings.Ts_st * diff(settings.index_T);
input.x0 = 0;
input.u = ones(1, settings.N);
input.x = [input.x0, input.x0 + cumsum(durations .* input.u)];
input.z = zeros(settings.nz, settings.N);
input.od = zeros(1, settings.N + 1);
input.y = zeros(1, settings.N);
input.yN = 0;
input.W = ones(1, settings.N);
input.WN = 0;
input.lambda = zeros(1, settings.N + 1);
input.mu_u = zeros(settings.N, 1);
input.lbu = -inf(1, settings.N);
input.ubu = inf(1, settings.N);
[input.lb, input.ub, input.lbx, input.ubx, input.mu_x, input.mu] = deal([]);
end

function tests = test_mex_workspace_lifecycle
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.assumeNotEmpty(mex.getCompilerConfigurations('C', 'Selected'), ...
    'Configure a MATLAB C compiler with mex -setup C to run the MEX lifecycle tests.');

testDir = fileparts(mfilename('fullpath'));
coreDir = fullfile(fileparts(testDir), 'mex_core');
sharedCore = fullfile(fileparts(testDir), 'core');
tmp = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(tmp.Folder));
testCase.addTeardown(@unloadFixtures);
options = {'-outdir', tmp.Folder, ['-I' coreDir]};

mex(options{:}, '-output', 'workspace_condensing', ...
    fullfile(coreDir, 'Condensing.c'), fullfile(coreDir, 'mpc_common.c'), '-lmwblas');
mex(options{:}, '-output', 'workspace_partial_condensing', ...
    fullfile(coreDir, 'partial_condensing_default.c'), ...
    fullfile(coreDir, 'partial_condensing_routines.c'), fullfile(coreDir, 'mpc_common.c'), '-lmwblas');
integratorSources = fullfile(coreDir, ...
    {'Simulate_System.c', 'sim.c', 'erk.c', 'irk_ode.c', 'irk_dae.c', 'mpc_common.c'});
mex(options{:}, ['-I' sharedCore], '-output', 'workspace_simulate', integratorSources{:}, ...
    fullfile(sharedCore, 'erk_kernel.c'), ...
    fullfile(testDir, 'fixtures', 'scalar_ode.c'), '-lmwblas', '-lmwlapack');
end

function unloadFixtures
clear workspace_condensing workspace_partial_condensing workspace_simulate
end

function testFullCondensingChangesHorizon(testCase)
for horizon = [2, 5, 3, 5]
    [mem, settings] = scalarQP(horizon);
    workspace_condensing(mem, settings);
    dynamics = tril(ones(horizon));
    expected = eye(horizon) + dynamics' * dynamics;
    testCase.verifyEqual(mem.Hc, expected, 'AbsTol', 1e-12);
end
end

function testPartialCondensingChangesBlockSize(testCase)
[mem, settings] = scalarQP(6);
for blockCount = [3, 2, 6, 3]
    settings.N2 = blockCount;
    [~, ~, Ccx, Ap, Bp, ap, lxc, uxc, Ccg, lgc, ugc] = workspace_partial_condensing(mem, settings);
    testCase.verifyEqual(Ap, ones(1, blockCount), 'AbsTol', 1e-12);
    testCase.verifyEqual(Bp, ones(1, settings.N), 'AbsTol', 1e-12);
    testCase.verifyEqual(ap, zeros(1, blockCount), 'AbsTol', 1e-12);
    testCase.verifySize(Ccx, [0, blockCount + settings.N]);
    testCase.verifySize(Ccg, [0, blockCount + settings.N]);
    testCase.verifyEqual([lxc, uxc, lgc, ugc], zeros(0, 4));
end
end

function testSimulatorChangesMethodStepAndStepCount(testCase)
[explicit, settings] = explicitConfig();

differentStep = explicit;
differentStep.h = 0.08;
differentSteps = differentStep;
differentSteps.num_steps = 3;

implicit = explicit;
implicit.sim_method = 2;
implicit.num_stages = 1;
implicit.A_tab = 0.5;
implicit.B_tab = 1;
implicit.newton_iter = 3;
implicit.JFK = implicit.h;

for configuration = {explicit, differentStep, differentSteps, implicit, explicit}
    config = configuration{1};
    [nextState, nextAlgebraic] = workspace_simulate(3, 2, [], 0, config, settings);
    testCase.verifyEqual(nextState, 3 + 2 * config.h * config.num_steps, 'AbsTol', 1e-12);
    testCase.verifySize(nextAlgebraic, [0, 1]);
end

dae = implicit;
dae.nz = 1;
dae.sim_method = 3;
dae.A_tab = 1;
dae.JFK = [dae.h, 0];
daeSettings = settings;
daeSettings.nz = 1;
[nextState, nextAlgebraic] = workspace_simulate(3, 2, 1, 0, dae, daeSettings);
testCase.verifyEqual(nextState, 3 + 2 * dae.h * dae.num_steps, 'AbsTol', 1e-12);
testCase.verifyEqual(nextAlgebraic, 0, 'AbsTol', 1e-12);
singleOutput = workspace_simulate(3, 2, 1, 0, dae, daeSettings);
testCase.verifyEqual(singleOutput, nextState, 'AbsTol', 1e-12);
[nextState, nextAlgebraic] = workspace_simulate(3, 2, [], 0, explicit, settings);
testCase.verifyEqual(nextState, 3 + 2 * explicit.h * explicit.num_steps, 'AbsTol', 1e-12);
testCase.verifySize(nextAlgebraic, [0, 1]);
end

function testSimulatorHandlesOneOrNoRequestedOutput(testCase)
[mem, settings] = explicitConfig();
expected = 3 + 2 * mem.h * mem.num_steps;
testCase.verifyEqual(workspace_simulate(3, 2, 0, [], mem, settings), expected, 'AbsTol', 1e-12);
testCase.verifyEqual(callWithNoRequestedOutput(3, 2, [], 0, mem, settings), expected, 'AbsTol', 1e-12);
end

function testSimulatorRejectsInvalidGatewayArguments(testCase)
[mem, settings] = explicitConfig();
args = {3, 2, [], 0, mem, settings};
argumentsError = 'MATMPC:SimulateSystem:Arguments';
inputError = 'MATMPC:SimulateSystem:Input';
configurationError = 'MATMPC:SimulateSystem:Configuration';
testCase.verifyError(@() workspace_simulate(args{1:5}), argumentsError);
testCase.verifyError(@() callWithThreeOutputs(args{:}), argumentsError);
for invalidState = {single(3), complex(3, 1), sparse(3), [3, 4]}
    testCase.verifyError(@() workspace_simulate(invalidState{1}, args{2:end}), inputError);
end
testCase.verifyError(@() workspace_simulate(3, [], [], 0, mem, settings), inputError);
testCase.verifyError(@() workspace_simulate(3, 2, [1, 2], 0, mem, settings), inputError);
testCase.verifyError(@() workspace_simulate(3, 2, [], 0, mem, rmfield(settings, 'nx')), configurationError);
settings.np = 2;
testCase.verifyError(@() workspace_simulate(3, 2, [], 0, mem, settings), inputError);

% A rejected call leaves the persistent integrator ready to use.
state = workspace_simulate(3, 2, [], [0; 0], mem, settings);
testCase.verifyEqual(state, 3 + 2 * mem.h * mem.num_steps, 'AbsTol', 1e-12);
end

function result = callWithNoRequestedOutput(varargin)
workspace_simulate(varargin{:});
result = ans; %#ok<NOANS>
end

function callWithThreeOutputs(varargin)
[~, ~, ~] = workspace_simulate(varargin{:});
end

function [mem, settings] = explicitConfig()
settings = struct('nx', 1, 'nu', 1, 'nz', 0, 'np', 0);
mem = settings;
mem.sim_method = 1;
mem.A_tab = [0, 0, 0, 0; 0.5, 0, 0, 0; 0, 0.5, 0, 0; 0, 0, 1, 0];
mem.B_tab = [1/6, 1/3, 1/3, 1/6];
mem.num_stages = 4;
mem.num_steps = 2;
mem.h = 0.05;
mem.Sx = 1;
mem.Su = 0;
end

function [mem, settings] = scalarQP(horizon)
settings = struct('nx', 1, 'nu', 1, 'nc', 0, 'ncN', 0, 'nbx', 0, 'N', horizon);
mem.Q = ones(1, horizon + 1);
mem.S = zeros(1, horizon);
mem.R = ones(1, horizon);
mem.A = ones(1, horizon);
mem.B = ones(1, horizon);
mem.gx = zeros(1, horizon + 1);
mem.gu = zeros(1, horizon);
mem.a = zeros(1, horizon);
mem.ds0 = 0;
mem.Hc = zeros(horizon);
mem.gc = zeros(horizon, 1);
mem.iter = 1;
mem.hot_start = 0;
emptyFields = {'Cx', 'Cgx', 'Cgu', 'CgN', 'lc', 'uc', 'lb_dx', 'ub_dx', ...
    'Ccg', 'Ccx', 'lcc', 'ucc', 'lxc', 'uxc', 'idxc'};
for field = emptyFields
    mem.(field{1}) = [];
end
end

function tests = test_matmpc_controller
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.addTeardown(@path, path);
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'mex_core'));
gateway = which('qp_generation');
assumeFalse(testCase, isempty(gateway), 'Build InvertedPendulum before running controller integration tests.');
buildRoot = fileparts(fileparts(gateway));
metadata = fullfile(buildRoot, 'data', 'settings.mat');
assumeTrue(testCase, exist(metadata, 'file') == 2, 'Generated model settings required.');
loaded = load(metadata, 'settings');
assumeEqual(testCase, loaded.settings.model, 'InvertedPendulum');
settings = loaded.settings;
settings.N = 40;
testCase.TestData.settings = settings;
if ispc
    platform = 'win64';
elseif ismac
    platform = 'mac';
else
    platform = 'linux';
end
solver = fullfile(buildRoot, 'solver', platform);
testCase.TestData.qpOASES = exist(fullfile(solver, 'qpoases', ['qpOASES_sequence.', mexext]), 'file') ~= 0 ...
    || exist('qpOASES_sequence', 'file') == 3;
testCase.TestData.OSQP = exist(fullfile(solver, 'osqp', ['osqp_mex.', mexext]), 'file') ~= 0 ...
    || exist('osqp_mex', 'file') == 3;
end

function testInterleavedQpOasesAndReset(testCase)
assumeTrue(testCase, testCase.TestData.qpOASES);
settings = testCase.TestData.settings;
settingsB = settings;
settingsB.N = 30;
measuredA = [0, 0.01, 0.02; pi, pi-0.05, pi-0.10; 0, 0, 0; 0, 0, 0];
measuredB = [-0.1, -0.11, -0.12; pi-0.10, pi-0.15, pi-0.20; 0, 0, 0; 0, 0, 0];
expectedA = runSequence(settings, measuredA);
expectedB = runSequence(settingsB, measuredB);
a = MATMPCController(settings);
b = MATMPCController(settingsB);
release = onCleanup(@() delete([a, b]));
for k = 1:3
    verifyEqual(testCase, a.step(measuredA(:, k)), expectedA(k), 'AbsTol', 1e-8);
    verifyEqual(testCase, b.step(measuredB(:, k)), expectedB(k), 'AbsTol', 1e-8);
end
a.reset();
verifyEqual(testCase, a.snapshot().sample, 0);
verifyEqual(testCase, a.step(measuredA(:, 1)), expectedA(1), 'AbsTol', 1e-8);
delete(b);
delete(b);
verifyEqual(testCase, a.step(measuredA(:, 2)), expectedA(2), 'AbsTol', 1e-8);
end

function testOptionalDebugSnapshotsOwnTheirArrays(testCase)
assumeTrue(testCase, testCase.TestData.qpOASES);
controller = MATMPCController(testCase.TestData.settings);
release = onCleanup(@() delete(controller));
x0 = [0; pi; 0; 0];
[~, first] = controller.step(x0);
snapshot = controller.snapshot(true);
priorGradient = snapshot.qp.gx + zeros(size(snapshot.qp.gx));
priorControl = first.u + zeros(size(first.u));
verifyFalse(testCase, isfield(controller.snapshot(), 'qp'));
reference = repmat([0.4; 0; 0; 0; 0], 1, snapshot.settings.N);
controller.step(x0, struct('y', reference, 'yN', [0.4; 0; 0; 0]));
current = controller.snapshot(true);
verifyEqual(testCase, current.input.y, reference);
verifyGreaterThan(testCase, norm(current.qp.gx - priorGradient), 0.1);
verifyEqual(testCase, snapshot.qp.gx, priorGradient);
verifyEqual(testCase, first.u, priorControl);
snapshot.input.x(:) = 99;
verifyNotEqual(testCase, controller.snapshot().input.x, snapshot.input.x);
end

function testFailureKeepsTheLastCommittedIterate(testCase)
assumeTrue(testCase, testCase.TestData.qpOASES);
settings = testCase.TestData.settings;
c = MATMPCController(settings);
other = MATMPCController(settings);
release = onCleanup(@() delete([c, other]));
x0 = [0; pi; 0; 0];
c.step(x0);
other.step(x0);
before = c.snapshot();
verifyError(testCase, @() c.step(x0, struct('N', 5)), 'MATMPC:UnsupportedControllerUpdate');
% The binary throws directly on hotstart infeasibility on some platforms.
failed = false;
try
    c.step([100; pi; 0; 0]);
catch problem
    failed = true;
    verifyTrue(testCase, strcmp(problem.identifier, 'MATMPC:QPSolverFailed') || ...
        ~isempty(strfind(problem.message, 'qpOASES')));
end
assertTrue(testCase, failed);
verifyEqual(testCase, c.snapshot(), before);
verifyTrue(testCase, all(isfinite(other.step(x0))));
verifyTrue(testCase, all(isfinite(c.step(x0))));
verifyEqual(testCase, c.snapshot().sample, before.sample + 1);
end

function testMoveBlockingSnapshotContainsStageAndReducedQP(testCase)
assumeTrue(testCase, testCase.TestData.qpOASES);
settings = testCase.TestData.settings;
settings.r = 3;
settings.index_T = [0, 1, 10, settings.N];
controller = MATMPCController(settings, struct('qpsolver', 'qpoases_mb'));
release = onCleanup(@() delete(controller)); %#ok<NASGU>
x0 = [0; pi; 0; 0];
firstControl = controller.step(x0);
snapshot = controller.snapshot(true);
verifySize(testCase, snapshot.qp.R, [settings.nu, settings.N*settings.nu]);
verifySize(testCase, snapshot.qp.R_block, [settings.nu, settings.r*settings.nu]);
verifySize(testCase, snapshot.qp.gu_block, [settings.nu, settings.r]);
verifySize(testCase, snapshot.qp.Hc_r, [settings.r*settings.nu, settings.r*settings.nu]);
verifyFalse(testCase, isfield(snapshot.qp, 'Hc'));
controller.reset();
verifyEqual(testCase, controller.step(x0), firstControl, 'AbsTol', 1e-8);
end

function testOSQPOwnsBothWorkspaceForms(testCase)
assumeTrue(testCase, testCase.TestData.OSQP);
settings = testCase.TestData.settings;
x0 = [0; pi; 0; 0];
for solver = {'osqp_sparse', 'osqp_partial_sparse'}
    a = MATMPCController(settings, struct('qpsolver', solver{1}));
    b = MATMPCController(settings, struct('qpsolver', solver{1}));
    release = onCleanup(@() delete([a, b]));
    first = a.step(x0);
    verifyEqual(testCase, b.step(x0), first, 'AbsTol', 1e-8);
    delete(b);
    a.reset();
    verifyEqual(testCase, a.step(x0), first, 'AbsTol', 1e-8);
    clear release;
end
end

function testBuildGuardCountsAllLiveControllers(testCase)
assumeTrue(testCase, testCase.TestData.qpOASES);
settings = testCase.TestData.settings;
a = MATMPCController(settings);
b = MATMPCController(settings);
release = onCleanup(@() delete([a, b]));
verifyError(testCase, @() MATMPCController.assertBuildAllowed(), 'MATMPC:ActiveControllers');
verifyError(testCase, @() Build_Model('InvertedPendulum'), 'MATMPC:ActiveControllers');
verifyError(testCase, @() Compile_Mex(), 'MATMPC:ActiveControllers');
root = fileparts(fileparts(mfilename('fullpath')));
verifyError(testCase, @() run(fullfile(root, 'Simulation.m')), 'MATMPC:ActiveControllers');
verifyError(testCase, @() run(fullfile(root, 'Initialization_Simulink.m')), 'MATMPC:ActiveControllers');
delete(a);
b.reset();
verifyError(testCase, @() MATMPCController.assertBuildAllowed(), 'MATMPC:ActiveControllers');
delete(b);
MATMPCController.assertBuildAllowed();
end

function controls = runSequence(settings, measurements)
controller = MATMPCController(settings);
release = onCleanup(@() delete(controller));
controls = zeros(1, size(measurements, 2));
for k = 1:numel(controls)
    controls(k) = controller.step(measurements(:, k));
end
end

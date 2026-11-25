function tests = test_qpoases_failure_handling
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testDir = fileparts(mfilename('fullpath'));
solverDir = fullfile(fileparts(testDir), 'nmpc');
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(solverDir));

stubs = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
qpStub = fullfile(stubs.Folder, 'qpOASES_sequence.m');
recoverStub = fullfile(stubs.Folder, 'Recover.m');
copyfile(fullfile(testDir, 'fixtures', 'qpOASES_sequence.txt'), qpStub);
copyfile(fullfile(testDir, 'fixtures', 'Recover.txt'), recoverStub);
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(stubs.Folder));
testCase.assertEqual(which('qpOASES_sequence'), qpStub);
testCase.assertEqual(which('Recover'), recoverStub);
end

function testDefaultRejectsInfeasibleInitialization(testCase)
[sizes, mem] = qpInput();
mem.qpoases_opt.exitflag = -2;
verifyRejected(testCase, @() mpc_qp_solve_qpoases(sizes, mem), -2, 42);
end

function testDefaultRejectsIterationLimitWithExistingHandle(testCase)
[sizes, mem] = qpInput();
mem.warm_start = 23;
mem.qpoases_opt.exitflag = 1;
verifyRejected(testCase, @() mpc_qp_solve_qpoases(sizes, mem), 1, []);
end

function testMoveBlockingRejectsInfeasibleInitialization(testCase)
[sizes, mem] = moveBlockingInput();
mem.qpoases_opt.exitflag = -2;
verifyRejected(testCase, @() mpc_qp_solve_qpoases_mb(sizes, mem, struct()), -2, 42);
end

function testMoveBlockingRejectsIterationLimitOnHotstart(testCase)
[sizes, mem] = moveBlockingInput();
mem.warm_start = 23;
mem.hot_start = 1;
mem.qpoases_opt.exitflag = 1;
verifyRejected(testCase, @() mpc_qp_solve_qpoases_mb(sizes, mem, struct()), 1, []);
end

function testDefaultSuccessStillRecovers(testCase)
[sizes, mem] = qpInput();
recoveries = Recover();
[time, mem] = mpc_qp_solve_qpoases(sizes, mem);
testCase.verifyEqual(Recover(), recoveries + 1);
testCase.verifyEqual(mem.warm_start, 42);
testCase.verifyEqual(mem.du, [0.25, -0.5]);
testCase.verifyEqual(mem.mu_u_new, [-1; 2]);
testCase.verifyEqual(time, 2);
end

function testMoveBlockingSuccessStillRecovers(testCase)
[sizes, mem] = moveBlockingInput();
recoveries = Recover();
[time, mem] = mpc_qp_solve_qpoases_mb(sizes, mem, struct());
testCase.verifyEqual(Recover(), recoveries + 1);
testCase.verifyEqual(mem.warm_start, 42);
testCase.verifyEqual(mem.du, [0.25, 0.25]);
testCase.verifyEqual(mem.mu_u_new, [-0.5; -0.5]);
testCase.verifyEqual(time, 2);
end

function verifyRejected(testCase, solve, exitflag, expectedClosedHandle)
recoveries = Recover();
closedHandles = qpOASES_sequence('testClosedHandles');

try
    solve();
catch failure
    testCase.verifyEqual(failure.identifier, 'MATMPC:QPSolverFailed');
    expectedText = sprintf('exitflag %d', exitflag);
    testCase.verifyThat(failure.message, matlab.unittest.constraints.ContainsSubstring(expectedText));
    testCase.verifyEqual(Recover(), recoveries);
    testCase.verifyEqual(qpOASES_sequence('testClosedHandles'), [closedHandles; expectedClosedHandle(:)]);
    return;
end
testCase.verifyFail('Expected the QP solve to be rejected.');
end

function [sizes, mem] = qpInput()
sizes = struct('N', 2, 'nu', 1, 'nbx', 0, 'nc', 0, 'ncN', 0);
mem.warm_start = 0;
mem.hot_start = 0;
mem.Hc = eye(2);
mem.gc = zeros(2, 1);
mem.Ccx = zeros(0, 2);
mem.Ccg = zeros(0, 2);
mem.lb_du = -ones(2, 1);
mem.ub_du = ones(2, 1);
mem.lxc = zeros(0, 1);
mem.uxc = zeros(0, 1);
mem.lcc = zeros(0, 1);
mem.ucc = zeros(0, 1);
mem.du = zeros(1, 2);
mem.mu_u_new = zeros(2, 1);
mem.mu_x_new = zeros(0, 1);
mem.mu_new = zeros(0, 1);
mem.qpoases_opt = struct('exitflag', 0, 'solution', [0.25; -0.5], 'multipliers', [1; -2]);
end

function [sizes, mem] = moveBlockingInput()
[sizes, mem] = qpInput();
mem.r = 1;
mem.index_T = [0, 2];
mem.T = [1; 1];
mem.Hc_r = 1;
mem.gc_r = 0;
mem.Ccx_r = zeros(0, 1);
mem.Ccg_r = zeros(0, 1);
mem.lb_du = [0.25; 0.25];
mem.qpoases_opt.solution = 0.25;
mem.qpoases_opt.multipliers = 1;
end

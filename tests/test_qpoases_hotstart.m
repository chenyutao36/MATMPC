function tests = test_qpoases_hotstart
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testDir = fileparts(mfilename('fullpath'));
root = fileparts(testDir);
testCase.assumeNotEmpty(mex.getCompilerConfigurations('C', 'Selected'), ...
    'Configure a MATLAB C compiler for the Recover and alias-test MEX files.');
if ismac
    platform = 'mac';
elseif isunix
    platform = 'linux';
else
    platform = 'win64';
end
solverDir = fullfile(root, 'solver', platform, 'qpoases');
native = fullfile(solverDir, ['qpOASES_sequence.', mexext]);
testCase.assumeTrue(exist(native, 'file') ~= 0, 'A native qpOASES sequence MEX file is required.');

testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(root, 'nmpc')));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(solverDir));
temporary = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
copyfile(native, fullfile(temporary.Folder, ['qpOASES_sequence_trace_native.', mexext]));
copyfile(fullfile(testDir, 'fixtures', 'qpOASES_sequence_trace.txt'), ...
    fullfile(temporary.Folder, 'qpOASES_sequence.m'));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(temporary.Folder));
testCase.addTeardown(@unloadFixtures);
testCase.assertEqual(exist('qpOASES_sequence_trace_native', 'file'), 3);
core = fullfile(root, 'core');
mex('-silent', '-outdir', temporary.Folder, ['-I', core], ...
    fullfile(root, 'mex_core', 'Recover.c'), fullfile(core, 'recover_kernel.c'));
mex('-silent', '-outdir', temporary.Folder, fullfile(testDir, 'fixtures', 'overwrite_qp_matrix_test.c'));
end

function testStandardOptInOnlyUsesShortCallForExactMatches(testCase)
exerciseUpdates(testCase, false, true, 'ihmmh');
end

function testMoveBlockedOptInOnlyUsesShortCallForExactMatches(testCase)
exerciseUpdates(testCase, true, true, 'ihmmh');
end

function testStandardDefaultsToMatrixUpdates(testCase)
exerciseUpdates(testCase, false, false, 'immmm');
end

function testMoveBlockedDefaultsToMatrixUpdates(testCase)
exerciseUpdates(testCase, true, false, 'immmm');
end

function exerciseUpdates(testCase, blocked, hotstart, expectedCalls)
[mem, solve, H, A, g, lift] = qpCase(blocked, hotstart);
qpOASES_sequence('testResetTrace');
[~, mem] = solve(mem);
sequence = mem.warm_start;
closeSequence = onCleanup(@() qpOASES_sequence('c', sequence));
testCase.verifyEqual(mem.du, lift([2; 1]), 'AbsTol', 1e-9);

mem.(g) = [-6; -8];
[~, mem] = solve(mem);
testCase.verifyEqual(mem.du, lift([2; 2]), 'AbsTol', 1e-9);

% Exercise the write-through pattern of legacy MEX code. Ordinary MATLAB
% assignments would detach data and could hide a shared-cache bug.
overwrite_qp_matrix_test(mem.(H), 4, 8);
testCase.verifyEqual(mem.(H)(2, 2), 8);
[~, mem] = solve(mem);
testCase.verifyEqual(mem.du, lift([2; 1]), 'AbsTol', 1e-9);

overwrite_qp_matrix_test(mem.(A), 1, 2);
testCase.verifyEqual(mem.(A)(1, 1), 2);
[~, mem] = solve(mem);
testCase.verifyEqual(mem.du, lift([1; 1]), 'AbsTol', 1e-9);
[~, mem] = solve(mem);
testCase.verifyEqual(mem.du, lift([1; 1]), 'AbsTol', 1e-9);
testCase.verifyEqual(qpOASES_sequence('testCommands'), expectedCalls);
testCase.verifyEqual(isfield(mem, 'qpoases_matrix_cache'), hotstart);
end

function testFailedMatrixUpdateInvalidatesPreviousValues(testCase)
[mem, solve, H, ~, ~, lift] = qpCase(false, true);
qpOASES_sequence('testResetTrace');
[~, mem] = solve(mem);
sequence = mem.warm_start;
closeSequence = onCleanup(@() qpOASES_sequence('c', sequence));
failed = mem;
failed.(H) = diag([8, 8]);
qpOASES_sequence('testFailAfterUpdate');
testCase.verifyError(@() solve(failed), 'MATMPC:QPSolverFailed');

% The fault occurred after the actual native sequence installed failed.Hc.
% Retrying the previous value must reinstall its original matrices.
[~, mem] = solve(mem);
testCase.verifyEqual(mem.du, lift([2; 1]), 'AbsTol', 1e-9);
testCase.verifyEqual(qpOASES_sequence('testCommands'), 'imm');
end

function testFlagChangesDoNotArmOlderAliases(testCase)
[mem, solve] = qpCase(false, false);
qpOASES_sequence('testResetTrace');
[~, mem] = solve(mem);
sequence = mem.warm_start;
closeSequence = onCleanup(@() qpOASES_sequence('c', sequence));
mem.hot_start = 1;
[~, mem] = solve(mem);
[~, mem] = solve(mem);
testCase.verifyFalse(isfield(mem, 'qpoases_matrix_cache'));
testCase.verifyEqual(qpOASES_sequence('testCommands'), 'imm');
clear closeSequence;

[mem, solve] = qpCase(false, true);
qpOASES_sequence('testResetTrace');
[~, mem] = solve(mem);
sequence = mem.warm_start;
closeSequence = onCleanup(@() qpOASES_sequence('c', sequence));
previous = mem;
mem.hot_start = 0;
[~, mem] = solve(mem);
[~, previous] = solve(previous);
[~, previous] = solve(previous);
testCase.verifyFalse(isfield(mem, 'qpoases_matrix_cache'));
testCase.verifyFalse(isfield(previous, 'qpoases_matrix_cache'));
testCase.verifyEqual(qpOASES_sequence('testCommands'), 'immm');
end

function [mem, solve, H, A, g, lift] = qpCase(blocked, hotstart)
N = 2;
H = 'Hc';
A = 'Ccx';
g = 'gc';
map = eye(2);
if blocked
    N = 4;
    H = 'Hc_r';
    A = 'Ccx_r';
    g = 'gc_r';
    map = [1 0; 1 0; 0 1; 0 1];
end
sizes = struct('nx', 1, 'nu', 1, 'N', N, 'nbx', 1, 'nc', 0, 'ncN', 0);
if blocked
    solve = @(state) mpc_qp_solve_qpoases_mb(sizes, state, struct());
    mem.r = 2;
    mem.index_T = [0 2 4];
    mem.T = map;
    mem.Ccg_r = zeros(0, 2);
else
    solve = @(state) mpc_qp_solve_qpoases(sizes, state);
    mem.Ccg = zeros(0, 2);
end
lift = @(u) reshape(map * u, 1, N);
mem.warm_start = 0;
mem.hot_start = double(hotstart);
mem.qpoases_opt = qpOASES_options('MPC');
mem.(H) = diag([2, 4]);
mem.(g) = [-6; -4];
mem.(A) = zeros(N, 2);
mem.(A)(1, 1) = 1;
mem.lb_du = -10 * ones(N, 1);
mem.ub_du = 10 * ones(N, 1);
mem.lxc = -inf(N, 1);
mem.uxc = [2; inf(N-1, 1)];
mem.lcc = zeros(0, 1);
mem.ucc = zeros(0, 1);
mem.Q = zeros(1, N+1);
mem.S = zeros(1, N);
mem.A = ones(1, N);
mem.B = zeros(1, N);
mem.Cx = 1;
mem.Cgx = zeros(0, N);
mem.CgN = [];
mem.gx = zeros(1, N+1);
mem.a = zeros(1, N);
mem.ds0 = 0;
mem.du = zeros(1, N);
mem.mu_u_new = zeros(N, 1);
mem.mu_x_new = zeros(N, 1);
mem.mu_new = zeros(0, 1);
end

function unloadFixtures()
clear Recover overwrite_qp_matrix_test qpOASES_sequence qpOASES_sequence_trace_native;
end

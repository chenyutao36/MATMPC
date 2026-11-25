function tests = test_qpoases_move_blocking
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
projectRoot = fileparts(fileparts(mfilename('fullpath')));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(projectRoot, 'nmpc')));
if ismac
    platform = 'mac';
elseif isunix
    platform = 'linux';
else
    platform = 'win64';
end
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(...
    fullfile(projectRoot, 'solver', platform, 'qpoases')));
testCase.assumeEqual(exist('qpOASES_sequence', 'file'), 3, ...
    'This test requires a native qpOASES sequence MEX file.');
testCase.assumeTrue(~isempty(mex.getCompilerConfigurations('C', 'Selected')), ...
    'This test requires a selected MATLAB C compiler to build Recover.');

build = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
source = fullfile(projectRoot, 'mex_core', 'Recover.c');
core = fullfile(projectRoot, 'core');
mex('-silent', '-outdir', build.Folder, ['-I', core], source, fullfile(core, 'recover_kernel.c'));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(build.Folder));
testCase.addTeardown(@clearRecover);
end

function testMultiInputBoundsAndDuals(testCase)
[sizes, mem] = sampleQP();
[~, mem] = mpc_qp_solve_qpoases_mb(sizes, mem, struct());
closeSequence = onCleanup(@() qpOASES_sequence('c', mem.warm_start));
verifySolution(testCase, mem, [0, 3, 0, 0; -3, 0, 0, 0]);

% This call updates the existing QP, including the Hessian and Jacobian.
mem.gc_r(1) = -5;
[~, mem] = mpc_qp_solve_qpoases_mb(sizes, mem, struct());
verifySolution(testCase, mem, [0, 4, 0, 0; -3, 0, 0, 0]);

% Opting in after initialization keeps the matrix-update path for this
% sequence. The lower input bound now ties across both stages of the block.
mem.hot_start = 1;
mem.lb_du(4) = -2;
[~, mem] = mpc_qp_solve_qpoases_mb(sizes, mem, struct());
verifySolution(testCase, mem, [0, 4, 0, 0; -1.5, -1.5, 0, 0]);
end

function testDisjointStageBoundsAreInfeasible(testCase)
[sizes, mem] = sampleQP();
% Both stages are feasible individually, but the common increment cannot
% satisfy the first stage's lower bound and the second stage's upper bound.
mem.lb_du(1) = 1.5;
testCase.verifyError(@() mpc_qp_solve_qpoases_mb(sizes, mem, struct()), ...
    'MATMPC:QPSolverFailed');
end

function verifySolution(testCase, mem, stageInputDuals)
tol = 1e-8;
testCase.verifyEqual(mem.du, [1, 1, 2, 2; -2, -2, -3, -3], 'AbsTol', tol);
testCase.verifyEqual(mem.mu_u_new, stageInputDuals(:), 'AbsTol', tol);
testCase.verifyEqual(mem.mu_x_new, [0; 0; 2; 0], 'AbsTol', tol);
testCase.verifyEqual(mem.mu_new, [0; 0; 0; 0; -3], 'AbsTol', tol);
testCase.verifyEqual(mem.lambda_new(end), -3, 'AbsTol', tol);
end

function [sizes, mem] = sampleQP()
sizes = struct('N', 4, 'nx', 1, 'nu', 2, 'nbx', 1, 'nc', 1, 'ncN', 1);
mem.r = 2;
mem.index_T = [0, 2, 4];
mem.T = kron([1, 0; 1, 0; 0, 1; 0, 1], eye(2));
mem.warm_start = 0;
mem.hot_start = 0;
mem.qpoases_opt = qpOASES_options('MPC');
mem.Hc_r = eye(4);
mem.gc_r = [-4; 5; -4; 6];
mem.lb_du = reshape([-10, -10, -10, -10; -2, -4, -10, -10], [], 1);
mem.ub_du = reshape([2, 1, 10, 10; 10, 10, 10, 10], [], 1);
mem.Ccx_r = zeros(4, 4);
mem.Ccx_r(3, 3) = 1;
mem.lxc = -inf(4, 1);
mem.uxc = [inf; inf; 2; inf];
mem.Ccg_r = zeros(5, 4);
mem.Ccg_r(5, 4) = 1;
mem.lcc = [-inf; -inf; -inf; -inf; -3];
mem.ucc = inf(5, 1);

% Supply the full-space buffers consumed by Recover.
mem.Q = zeros(1, 5);
mem.S = zeros(1, 8);
mem.R = zeros(2, 8);
mem.A = ones(1, 4);
mem.B = zeros(1, 8);
mem.Cx = 1;
mem.Cgx = zeros(1, 4);
mem.Cgu = zeros(1, 8);
mem.CgN = 1;
mem.gx = zeros(1, 5);
mem.a = zeros(1, 4);
mem.ds0 = 0;
mem.dx = zeros(1, 5);
mem.du = zeros(2, 4);
mem.mu_u_new = zeros(8, 1);
mem.mu_x_new = zeros(4, 1);
mem.mu_new = zeros(5, 1);
mem.lambda_new = zeros(1, 5);
end

function clearRecover()
clear Recover;
end

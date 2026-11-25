function tests = test_recover_mex
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
root = fileparts(fileparts(mfilename('fullpath')));
testCase.assumeTrue(~isempty(mex.getCompilerConfigurations('C', 'Selected')), ...
    'The Recover gateway test requires a selected MATLAB C compiler.');
build = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
source = fullfile(root, 'mex_core', 'Recover.c');
core = fullfile(root, 'core');
mex('-silent', '-outdir', build.Folder, ['-I', core], source, fullfile(core, 'recover_kernel.c'));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(build.Folder));
testCase.addTeardown(@clearRecover);
end

function testNumericalParityAndAliasIsolation(testCase)
[mem, settings] = sampleInputs();
alias = mem;
[mem.dx, mem.lambda_new] = Recover(mem, settings);
verifyTrajectory(testCase, mem.dx, mem.lambda_new);
testCase.verifyEqual(alias.dx, 17 * ones(2, 3));
testCase.verifyEqual(alias.lambda_new, -31 * ones(2, 3));
testCase.verifyEqual(alias.du, [0.5 -2; -1 0.5; 2 1]);
end

function testAcceptsFlatControlAndNeedsNoExistingOutputFields(testCase)
[mem, settings] = sampleInputs();
mem = rmfield(mem, {'dx', 'lambda_new'});
mem.du = mem.du(:);
[dx, lambda] = Recover(mem, settings);
verifyTrajectory(testCase, dx, lambda);
end

function testRejectsWrongArityAndTypes(testCase)
[mem, settings] = sampleInputs();
testCase.verifyError(@() noOutputs(mem, settings), 'MATMPC:Recover:InvalidArgumentCount');
testCase.verifyError(@() oneOutput(mem, settings), 'MATMPC:Recover:InvalidArgumentCount');
testCase.verifyError(@() twoOutputs(mem), 'MATMPC:Recover:InvalidArgumentCount');
testCase.verifyError(@() twoOutputs([], settings), 'MATMPC:Recover:InvalidStructure');
testCase.verifyError(@() twoOutputs([mem mem], settings), 'MATMPC:Recover:InvalidStructure');
end

function testRejectsInvalidAndOverflowingDimensions(testCase)
[mem, settings] = sampleInputs();
testCase.verifyError(@() twoOutputs(mem, rmfield(settings, 'N')), ...
    'MATMPC:Recover:InvalidDimensions');
for horizon = [0, -1, 1.5, Inf, NaN, 2^60]
    settings.N = horizon;
    testCase.verifyError(@() twoOutputs(mem, settings), 'MATMPC:Recover:InvalidDimensions');
end
end

function testRejectsMissingOrInvalidBuffers(testCase)
[mem, settings] = sampleInputs();
testCase.verifyError(@() twoOutputs(rmfield(mem, 'Q'), settings), 'MATMPC:Recover:InvalidBuffer');
invalid = mem;
invalid.Q = invalid.Q(:);
testCase.verifyError(@() twoOutputs(invalid, settings), 'MATMPC:Recover:InvalidBuffer');
invalid = mem;
invalid.B = single(invalid.B);
testCase.verifyError(@() twoOutputs(invalid, settings), 'MATMPC:Recover:InvalidBuffer');
invalid = mem;
invalid.Cgx = sparse(invalid.Cgx);
testCase.verifyError(@() twoOutputs(invalid, settings), 'MATMPC:Recover:InvalidBuffer');
invalid = mem;
invalid.mu_new = invalid.mu_new + 1i;
testCase.verifyError(@() twoOutputs(invalid, settings), 'MATMPC:Recover:InvalidBuffer');
invalid = mem;
invalid.mu_x_new = invalid.mu_x_new(1);
testCase.verifyError(@() twoOutputs(invalid, settings), 'MATMPC:Recover:InvalidBuffer');
end

function verifyTrajectory(testCase, dx, lambda)
testCase.verifyEqual(dx, [1 2 1; -2 -2 -5], 'AbsTol', 1e-12);
testCase.verifyEqual(lambda, [11.375 -10.125 0.75; -47.5 -66 -25], 'AbsTol', 1e-12);
end

function [mem, settings] = sampleInputs()
settings = struct('nx', 2, 'nu', 3, 'N', 2, 'nbx', 1, 'nc', 2, 'ncN', 1);
mem.Q = [2 1 4 -1 6 2; 1 3 -1 2 2 5];
mem.S = [1 2 0 0 1 2; -1 0 3 2 -1 0];
mem.A = [1 2 0.5 0; 0 -1 1 2];
mem.B = [1 0 2 0 2 -2; -1 3 0 1 0 1];
mem.Cx = [2 -1];
mem.Cgx = [1 2 0 3; -1 1 2 -2];
mem.CgN = [1 -2];
mem.gx = [1 -2 4; -1 3 1];
mem.a = [0.5 1; -0.5 -2];
mem.ds0 = [1; -2];
mem.du = [0.5 -2; -1 0.5; 2 1];
mem.mu_new = [0.5; 1.25; -2; 0.5; 1.75];
mem.mu_x_new = [1.5; -0.5];
mem.dx = 17 * ones(2, 3);
mem.lambda_new = -31 * ones(2, 3);
end

function noOutputs(varargin)
Recover(varargin{:});
end

function oneOutput(varargin)
result = Recover(varargin{:});
end

function twoOutputs(varargin)
[~, ~] = Recover(varargin{:});
end

function clearRecover()
clear Recover;
end

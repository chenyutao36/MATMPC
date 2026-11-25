function tests = test_condense_mb_mex
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.assumeNotEmpty(mex.getCompilerConfigurations('C', 'Selected'), ...
    'Configure a MATLAB C compiler to test the move-blocking gateway.');
projectRoot = fileparts(fileparts(mfilename('fullpath')));
core = fullfile(projectRoot, 'core');
gateway = fullfile(projectRoot, 'mex_core');
temporary = testCase.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(temporary.Folder));
testCase.applyFixture(matlab.unittest.fixtures.PathFixture(fullfile(projectRoot, 'nmpc')));
testCase.addTeardown(@() clearCondensing());
mex('-silent', '-outdir', temporary.Folder, ['-I' core], ...
    fullfile(gateway, 'Condensing_mb.c'), fullfile(core, 'condense_mb.c'), ...
    fullfile(core, 'condense_affine.c'), fullfile(gateway, 'matmpc_la_mex.c'), '-lmwblas');
end

function testReducedQPMatchesDirectSubstitution(testCase)
settings = struct('nx', 2, 'nu', 2, 'N', 4, 'nbx', 2, 'nc', 1, 'ncN', 2);
N = settings.N;
nx = settings.nx;
nu = settings.nu;
mem.Q = repmat([3 80; 1 2], 1, N+1);
mem.S = reshape(sin(1:16), nx, N*nu);
mem.R = repmat([2 0.5; 0.5 3], 1, N);
mem.A = reshape(cos(1:16), nx, N*nx);
mem.B = reshape(sin(2:17), nx, N*nu);
mem.Cx = [1 0.5; -0.25 2];
mem.Cgx = reshape(cos(2:9), 1, N*nx);
mem.Cgu = reshape(sin(4:11), 1, N*nu);
mem.CgN = [1 2; -1 3];
mem.gx = reshape(1:10, nx, N+1) / 10;
mem.gu = reshape(2:9, nu, N) / 10;
mem.a = reshape(3:10, nx, N) / 4;
mem.ds0 = [1; -2];
mem.lc = (1:6)' / 3;
mem.uc = 1 + mem.lc;
mem.lb_dx = -(1:8)' / 2;
mem.ub_dx = 1 + (1:8)' / 3;
mem.Hc_r = -17 * ones(4);
mem.r = 2;

for partition = {[0 1 4], [0 3 4], [0 2 4], [0 1 4]}
    mem.index_T = partition{1};
    mem.R_block = zeros(nu, mem.r*nu);
    mem.gu_block = zeros(nu, mem.r);
    T = zeros(N, mem.r);
    for block = 1:mem.r
        for stage = mem.index_T(block)+1:mem.index_T(block+1)
            reduced = (block-1)*nu + (1:nu);
            columns = (stage-1)*nu + (1:nu);
            mem.R_block(:, reduced) = mem.R_block(:, reduced) + mem.R(:, columns);
            mem.gu_block(:, block) = mem.gu_block(:, block) + mem.gu(:, stage);
            T(stage, block) = 1;
        end
    end
    expected = referenceQP(mem, settings, kron(T, eye(nu)));
    alias = mem;
    result = condense_move_blocked_qp(mem, settings);
    for field = fieldnames(expected)'
        name = field{1};
        testCase.verifyEqual(result.(name), expected.(name), 'AbsTol', 2e-12);
    end
    testCase.verifyEqual(result.Hc_r, result.Hc_r');
    testCase.verifyEqual(alias.Hc_r, -17 * ones(4));
end
end

function testMixedTermAndLegacyCall(testCase)
settings = struct('nx', 1, 'nu', 1, 'N', 2, 'nbx', 0, 'nc', 0, 'ncN', 0);
mem = struct('r', 1, 'index_T', [0 2], 'Q', [0 1 0], 'S', [0 1], ...
    'A', [0 0], 'B', [1 0], 'R_block', 1, 'gu_block', 0, 'gx', [0 0 0], ...
    'a', [0 0], 'ds0', 0, 'Cx', [], 'Cgx', [], 'Cgu', [], 'CgN', [], ...
    'lc', [], 'uc', [], 'lb_dx', [], 'ub_dx', [], ...
    'Hc_r', -8, 'gc_r', -13, 'Ccx_r', [], 'Ccg_r', [], ...
    'lcc', [], 'ucc', [], 'lxc', [], 'uxc', []);
qp = Condensing_mb(mem, settings);
testCase.verifyEqual(qp.Hc_r, 4);
testCase.verifyEqual(qp.gc_r, 0);
testCase.verifySize(qp.Ccx_r, [0 1]);
testCase.verifyEqual(mem.Hc_r, -8);

for invalid = {[0 0], [0 1], [0 NaN], [0 1.5]}
    wrong = mem;
    wrong.index_T = invalid{1};
    testCase.verifyError(@() Condensing_mb(wrong, settings), 'MATMPC:CondensingMB:Dimensions');
end
Condensing_mb(mem, settings);
testCase.verifyEqual(mem.Hc_r, 4);
testCase.verifyEqual(mem.gc_r, 0);
end

function expected = referenceQP(mem, settings, T)
n = size(T, 2);
eyeN = eye(n);
[base, generalOffset, stateOffset] = fullQP(mem, settings, zeros(size(T, 1), 1));
expected.Hc_r = zeros(n);
expected.gc_r = zeros(n, 1);
expected.Ccx_r = zeros(numel(stateOffset), n);
expected.Ccg_r = zeros(numel(generalOffset), n);
forward = zeros(n, 1);
for i = 1:n
    step = T * eyeN(:, i);
    [forward(i), general, state] = fullQP(mem, settings, step);
    backward = fullQP(mem, settings, -step);
    expected.gc_r(i) = (forward(i) - backward)/2;
    expected.Hc_r(i, i) = forward(i) + backward - 2*base;
    expected.Ccx_r(:, i) = state - stateOffset;
    expected.Ccg_r(:, i) = general - generalOffset;
    for j = 1:i-1
        both = fullQP(mem, settings, step + T * eyeN(:, j));
        expected.Hc_r(i, j) = both - forward(i) - forward(j) + base;
        expected.Hc_r(j, i) = expected.Hc_r(i, j);
    end
end
expected.lcc = mem.lc - generalOffset;
expected.ucc = mem.uc - generalOffset;
expected.lxc = mem.lb_dx - stateOffset;
expected.uxc = mem.ub_dx - stateOffset;
end

function [value, path, bounds] = fullQP(mem, sizes, controls)
x = mem.ds0;
value = 0;
path = zeros(sizes.N * sizes.nc + sizes.ncN, 1);
bounds = zeros(sizes.N * sizes.nbx, 1);
for i = 1:sizes.N
    xCols = (i-1)*sizes.nx + (1:sizes.nx);
    uCols = (i-1)*sizes.nu + (1:sizes.nu);
    u = controls(uCols);
    Q = tril(mem.Q(:, xCols));
    Q = Q + tril(Q, -1)';
    value = value + 0.5*(x'*Q*x + u'*mem.R(:, uCols)*u) + ...
        x'*mem.S(:, uCols)*u + mem.gx(:, i)'*x + mem.gu(:, i)'*u;
    path((i-1)*sizes.nc + (1:sizes.nc)) = mem.Cgx(:, xCols)*x + mem.Cgu(:, uCols)*u;
    x = mem.a(:, i) + mem.A(:, xCols)*x + mem.B(:, uCols)*u;
    bounds((i-1)*sizes.nbx + (1:sizes.nbx)) = mem.Cx*x;
end
Q = tril(mem.Q(:, sizes.N*sizes.nx + (1:sizes.nx)));
Q = Q + tril(Q, -1)';
value = value + 0.5*x'*Q*x + mem.gx(:, end)'*x;
path(sizes.N*sizes.nc+1:end) = mem.CgN*x;
end

function clearCondensing()
clear Condensing_mb;
end

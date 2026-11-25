function tests = test_sparse_qp_assembly
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    testCase.addTeardown(@path, path);
    root = fileparts(fileparts(mfilename('fullpath')));
    addpath(fullfile(root, 'nmpc'));
    if ispc
        platform = 'win64';
    elseif ismac
        platform = 'mac';
    else
        platform = 'linux';
    end
    solverDir = fullfile(root, 'solver', platform, 'osqp');
    available = exist(fullfile(solverDir, ['osqp_mex.', mexext]), 'file') ~= 0;
    if available
        addpath(solverDir);
    end
    testCase.TestData.hasOSQP = available;
end

function testAssemblyMatchesDenseReference(testCase)
    cases = {settings(4, 2, 0, 0, 0), settings(4, 2, 2, 1, [3, 1]), ...
        settings(2, 4, 6, 1, [3, 1])};
    for k = 1:numel(cases)
        sz = cases{k};
        mem = coefficients(sz);
        [H, D, q, l, u, G, C] = denseReference(mem, sz);
        [mem.sparse_qp, pPattern, aPattern] = sparse_qp_layout(sz, true);
        [p, a, qActual, lActual, uActual] = sparse_qp_data(mem);
        verifyTrue(testCase, issparse(pPattern) && issparse(aPattern));
        verifyEqual(testCase, p, H(find(pPattern)));
        verifyEqual(testCase, a, D(find(aPattern)));
        verifyEqual(testCase, qActual, q);
        verifyEqual(testCase, lActual, l);
        verifyEqual(testCase, uActual, u);

        mem.sparse_qp = sparse_qp_layout(sz, false);
        [hActual, aActual, qActual, lActual, uActual] = sparse_qp_data(mem);
        verifyEqual(testCase, hActual, sparse(H));
        verifyEqual(testCase, aActual, sparse([C; G]));
        verifyEqual(testCase, qActual, q);
        verifyEqual(testCase, lActual, [mem.lc; mem.ds0; -mem.a(:)]);
        verifyEqual(testCase, uActual, [mem.uc; mem.ds0; -mem.a(:)]);
    end
end

function testRealOSQPUpdatesPreserveStructuralZeros(testCase)
    assumeTrue(testCase, testCase.TestData.hasOSQP, 'Bundled OSQP MEX unavailable.');
    sz = settings(5, 2, 2, 1, [3, 1]);
    mem = coefficients(sz);
    mem.Q = repmat(2 * eye(sz.nx), 1, sz.N + 1);
    mem.R = repmat(2 * eye(sz.nu), 1, sz.N);
    mem.S(:) = 0;
    mem.Cgx(:) = 0;
    [mem.sparse_qp, pPattern, aPattern] = sparse_qp_layout(sz, true);
    solver = osqp;
    options = osqp.default_settings();
    options.eps_abs = 1e-7;
    options.eps_rel = 1e-7;
    options.verbose = false;
    limits = inf(size(aPattern, 1), 1);
    solver.setup(pPattern, [], aPattern, -limits, limits, options);

    for iteration = 1:2
        if iteration == 2
            mem.S(1) = 0.3;
            mem.Cgx(1) = 0.2;
        end
        [p, a, q, l, u] = sparse_qp_data(mem);
        solver.update('Px', p, 'Ax', a, 'q', q, 'l', l, 'u', u);
        actual = solver.solve();
        [H, D, q, l, u] = denseReference(mem, sz);
        reference = osqp;
        reference.setup(sparse(H), q, sparse(D), l, u, options);
        expected = reference.solve();
        assertEqual(testCase, actual.info.status, 'solved');
        assertEqual(testCase, expected.info.status, 'solved');
        verifyEqual(testCase, actual.x, expected.x, 'AbsTol', 2e-6);
        verifyEqual(testCase, actual.y, expected.y, 'AbsTol', 2e-6);
    end
end

function testStructureScalesWithHorizon(testCase)
    base = settings(250, 2, 2, 1, [3, 1]);
    [layout, p, a] = sparse_qp_layout(base, true);
    small = whos('layout', 'p', 'a');
    base.N = 1000;
    [layout, p, a] = sparse_qp_layout(base, true);
    large = whos('layout', 'p', 'a');
    verifyLessThan(testCase, sum([large.bytes]), 4.1 * sum([small.bytes]));
    denseHessianBytes = 8 * ((base.N + 1) * base.nx + base.N * base.nu)^2;
    verifyLessThan(testCase, sum([large.bytes]), 0.01 * denseHessianBytes);
end

function sz = settings(N, nu, nc, ncN, indices)
    sz = struct('N', N, 'nx', 3, 'nu', nu, 'nc', nc, 'ncN', ncN, ...
        'nbx', nnz(indices), 'nbx_idx', indices);
end

function m = coefficients(s)
    m.Q = fixture(s.nx, s.nx*(s.N+1), 2);
    m.S = fixture(s.nx, s.nu*s.N, 3);
    m.R = fixture(s.nu, s.nu*s.N, 4);
    m.A = fixture(s.nx, s.nx*s.N, 5);
    m.B = fixture(s.nx, s.nu*s.N, 6);
    m.Cgx = fixture(s.nc, s.nx*s.N, 7);
    m.Cgu = fixture(s.nc, s.nu*s.N, 8);
    m.CgN = fixture(s.ncN, s.nx, 9);
    m.gx = fixture(s.nx, s.N+1, 10);
    m.gu = fixture(s.nu, s.N, 11);
    m.ds0 = fixture(s.nx, 1, 12);
    m.a = fixture(s.nx, s.N, 13);
    m.lb_du = -ones(s.N*s.nu, 1);
    m.ub_du = 2*ones(s.N*s.nu, 1);
    m.lb_dx = -3*ones(s.N*s.nbx, 1);
    m.ub_dx = 4*ones(s.N*s.nbx, 1);
    m.lc = -5*ones(s.N*s.nc+s.ncN, 1);
    m.uc = 6*ones(s.N*s.nc+s.ncN, 1);
end

function M = fixture(n, m, seed)
    M = reshape((mod((1:n*m) * seed, 17) - 8) / 50, n, m);
end

function [H, D, q, l, u, G, C] = denseReference(mem, s)
    neq = (s.N+1)*s.nx;
    nw = neq + s.N*s.nu;
    H = zeros(nw);
    G = zeros(neq, nw);
    C = zeros(s.N*s.nc+s.ncN, nw);
    U = zeros(s.N*s.nu, nw);
    X = zeros(s.N*s.nbx, nw);
    G(1:s.nx, 1:s.nx) = eye(s.nx);
    for k = 1:s.N
        x = (k-1)*s.nx + (1:s.nx);
        xNext = x+s.nx;
        u = (k-1)*s.nu + (1:s.nu);
        j = neq+u;
        c = (k-1)*s.nc + (1:s.nc);
        bx = (k-1)*s.nbx + (1:s.nbx);
        H(x, x) = mem.Q(:, x);
        H(x, j) = mem.S(:, u);
        H(j, x) = mem.S(:, u)';
        H(j, j) = mem.R(:, u);
        G(xNext, x) = mem.A(:, x);
        G(xNext, xNext) = -eye(s.nx);
        G(xNext, j) = mem.B(:, u);
        U(u, j) = eye(s.nu);
        X(bx, k*s.nx + s.nbx_idx(1:s.nbx)) = eye(s.nbx);
        C(c, x) = mem.Cgx(:, x);
        C(c, j) = mem.Cgu(:, u);
    end
    xTerminal = s.N*s.nx + (1:s.nx);
    H(xTerminal, xTerminal) = mem.Q(:, xTerminal);
    C(s.N*s.nc + (1:s.ncN), xTerminal) = mem.CgN;
    D = [G; U; X; C];
    q = [mem.gx(:); mem.gu(:)];
    l = [mem.ds0; -mem.a(:); mem.lb_du; mem.lb_dx; mem.lc];
    u = [mem.ds0; -mem.a(:); mem.ub_du; mem.ub_dx; mem.uc];
end

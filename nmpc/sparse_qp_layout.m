function [layout, P, A] = sparse_qp_layout(settings, useOSQP)
%SPARSE_QP_LAYOUT Build sparse coordinates for the stagewise QP.
% The variables are stacked states followed by controls. For OSQP the rows
% are dynamics, input bounds, state bounds, then general constraints.

    N = settings.N;
    nx = settings.nx;
    nu = settings.nu;
    nc = settings.nc;
    nbx = settings.nbx;
    neq = (N + 1) * nx;
    nineq = N * nc + settings.ncN;
    nw = neq + N * nu;
    x = reshape(1:neq, nx, N + 1);
    u = reshape(neq + (1:N * nu), nu, N);
    cx = reshape(1:N * nc, nc, N);
    terminal = N * nc + (1:settings.ncN);

    [qRow, qCol] = blockCoordinates(x, x);
    [sRow, sCol] = blockCoordinates(x(:, 1:N), u);
    [rRow, rCol] = blockCoordinates(u, u);
    hRow = [qRow; sRow; sCol; rRow];
    hCol = [qCol; sCol; sRow; rCol];

    [gRow, gCol] = blockCoordinates(x(:, 2:end), x(:, 1:N));
    [bRow, bCol] = blockCoordinates(x(:, 2:end), u);
    [cxRow, cxCol] = blockCoordinates(cx, x(:, 1:N));
    [cuRow, cuCol] = blockCoordinates(cx, u);
    [cnRow, cnCol] = blockCoordinates(terminal(:), x(:, end));

    diagRows = (1:neq)';
    diagCols = diagRows;
    fixed = [ones(nx, 1); -ones(N * nx, 1)];
    if useOSQP
        controlRows = neq + (1:N * nu)';
        stateRows = neq + N * nu + (1:N * nbx)';
        boundedStates = x(settings.nbx_idx(1:nbx), 2:end);
        diagRows = [diagRows; controlRows; stateRows];
        diagCols = [diagCols; u(:); boundedStates(:)];
        fixed = [fixed; ones(N * (nu + nbx), 1)];
        inequalityOffset = neq + N * (nu + nbx);
        equalityOffset = 0;
    else
        % Ipopt uses variable bounds separately and puts inequalities first.
        inequalityOffset = 0;
        equalityOffset = nineq;
        diagRows = diagRows + equalityOffset;
    end

    aRow = [diagRows; gRow + equalityOffset; bRow + equalityOffset; ...
        cxRow + inequalityOffset; cuRow + inequalityOffset; cnRow + inequalityOffset];
    aCol = [diagCols; gCol; bCol; cxCol; cuCol; cnCol];
    layout.useOSQP = useOSQP;
    layout.nw = nw;
    layout.na = nineq + neq + double(useOSQP) * N * (nu + nbx);
    layout.fixed = fixed;

    if useOSQP
        % MATLAB sparse matrices and OSQP store nonzeros column by column.
        % Cache the update order from structural positions, including positions
        % whose coefficient is zero at a particular iteration.
        upper = find(hRow <= hCol);
        [~, pOrder] = sortrows([hCol(upper), hRow(upper)]);
        [~, layout.aOrder] = sortrows([aCol, aRow]);
        layout.pOrder = upper(pOrder);

        P = sparse(hRow(upper), hCol(upper), ones(numel(upper), 1), nw, nw);
        A = sparse(aRow, aCol, [fixed; ones(numel(aRow) - numel(fixed), 1)], layout.na, nw);
    else
        layout.hRow = hRow;
        layout.hCol = hCol;
        layout.aRow = aRow;
        layout.aCol = aCol;
    end
end

function [i, j] = blockCoordinates(rows, cols)
    nr = size(rows, 1);
    nc = size(cols, 1);
    n = size(rows, 2);
    i = repmat(reshape(rows, nr, 1, n), 1, nc, 1);
    j = repmat(reshape(cols, 1, nc, n), nr, 1, 1);
    i = i(:);
    j = j(:);
end

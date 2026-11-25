function [P, A, q, l, u] = sparse_qp_data(mem)
%SPARSE_QP_DATA Assemble a QP from the current coefficients and cached layout.
% OSQP receives fixed-length CSC value vectors. Ipopt receives sparse matrices.

    layout = mem.sparse_qp;
    h = [mem.Q(:); mem.S(:); mem.S(:); mem.R(:)];
    a = [layout.fixed; mem.A(:); mem.B(:); mem.Cgx(:); mem.Cgu(:); mem.CgN(:)];
    q = [mem.gx(:); mem.gu(:)];
    dynamics = [mem.ds0(:); -mem.a(:)];

    if layout.useOSQP
        P = h(layout.pOrder);
        A = a(layout.aOrder);
        l = [dynamics; mem.lb_du(:); mem.lb_dx(:); mem.lc(:)];
        u = [dynamics; mem.ub_du(:); mem.ub_dx(:); mem.uc(:)];
    else
        P = sparse(layout.hRow, layout.hCol, h, layout.nw, layout.nw);
        A = sparse(layout.aRow, layout.aCol, a, layout.na, layout.nw);
        l = [mem.lc(:); dynamics];
        u = [mem.uc(:); dynamics];
    end
end

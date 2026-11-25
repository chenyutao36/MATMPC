function [cpt_qp, mem] = mpc_qp_solve_osqp(sizes,mem)
    
    N=sizes.N;  
    nbx=sizes.nbx;
    nx=sizes.nx;
    nu=sizes.nu;
    neq = (N+1)*nx;
    nw = (N+1)*nx+N*nu;
    
    [P, A, q, l, u] = sparse_qp_data(mem);
    mem.qp_obj.update('Px', P, 'Ax', A, 'q', q, 'l', l, 'u', u);
    sol = mem.qp_obj.solve();
            
    assert(strcmp(sol.info.status,'solved'), ['QP Error: ' sol.info.status]);

    mem.dx(:) = reshape(sol.x(1:neq,1),[nx,N+1]);
    mem.du(:) = reshape(sol.x(neq+1:nw,1),[nu,N]);
        
    mem.lambda_new(:) = reshape(sol.y(1:neq,1),[nx,N+1]);
    mem.mu_u_new(:) = sol.y(neq+1:neq+N*nu,1);
    mem.mu_x_new(:) = sol.y(neq+N*nu+1:neq+N*nu+N*nbx,1);
    mem.mu_new(:) = sol.y(neq+N*nu+N*nbx+1:end,1);
    
    cpt_qp   = sol.info.run_time*1e3;
               
end

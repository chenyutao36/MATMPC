function [cpt_qp, mem] = mpc_qp_solve_qpalm_sparse(sizes,mem)
    
    N=sizes.N;  
    nbx=sizes.nbx;
    nx=sizes.nx;
    nu=sizes.nu;
    nc=sizes.nc;
    ncN=sizes.ncN;
    neq = (N+1)*nx;
    nw = (N+1)*nx+N*nu;
    nineq = N*nc+ncN;
    
    full2sparse(mem,sizes);    
    A = [mem.sparse_dG; mem.sparse_dBu; mem.sparse_dBx; mem.sparse_dB];
    l = [-mem.sparse_G; mem.sparse_lb];
    u = [-mem.sparse_G; mem.sparse_ub];
        
    mem.qpalm_solver=qpalm;
    mem.qpalm_settings=mem.qpalm_solver.default_settings();

    mem.qpalm_solver.setup(mem.sparse_H,mem.sparse_g, A, l, u, ...
            [reshape(mem.dx,[(N+1)*nx,1]);reshape(mem.du,[N*nu,1])], [reshape(mem.lambda_new,[neq,1]);mem.mu_u_new;mem.mu_x_new;mem.mu_new], mem.qpalm_settings); 
    S = 'mem.qpalm_solver.solve()';
    [T,sol] = evalc(S);
            
    assert(strcmp(sol.info.status,'solved'), ['QP Error: ' sol.info.status]);

    mem.dx(:) = reshape(sol.x(1:neq,1),[nx,N+1]);
    mem.du(:) = reshape(sol.x(neq+1:nw,1),[nu,N]);
        
    mem.lambda_new(:) = reshape(sol.y(1:neq,1),[nx,N+1]);
    mem.mu_u_new(:) = sol.y(neq+1:neq+N*nu,1);
    mem.mu_x_new(:) = sol.y(neq+N*nu+1:neq+N*nu+N*nbx,1);
    mem.mu_new(:) = sol.y(neq+N*nu+N*nbx+1:end,1);
    
    cpt_qp   = sol.info.run_time*1e3;
    
    mem.qpalm_solver.delete();
               
end

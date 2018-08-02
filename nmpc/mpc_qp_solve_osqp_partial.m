function [cpt_qp, mem] = mpc_qp_solve_osqp_partial(settings,settings2,mem, mem2)
    
    N=settings2.N;  
    Nc=settings2.Nc;
    nbx=settings2.nbx;
    nbx_idx=settings2.nbx_idx;
    nx=settings2.nx;
    nu=settings2.nu;
    nc=settings2.nc;
    ncN=settings2.ncN;
    neq = (N+1)*nx;
    nw = (N+1)*nx+N*nu;
    
    full2sparse(mem2,settings2);    
    A = [mem2.sparse_dG; mem2.sparse_dBu; mem2.sparse_dBx; mem2.sparse_dB];
    l = [-mem2.sparse_G; mem2.sparse_lb];
    u = [-mem2.sparse_G; mem2.sparse_ub];
        
    v_H = mem.mem2.sparse_H(mem.mem2.sparse_H_idx);
    v_A = A(mem.mem2.sparse_A_idx);
    mem.mem2.qp_obj.update('Px',v_H,'Ax',v_A,'q',mem.mem2.sparse_g,'l',l,'u',u);    
    sol = mem.mem2.qp_obj.solve();
    
    assert(strcmp(sol.info.status,'solved'), ['QP Error: QP is ' sol.info.status]);

    mem.du(:) = reshape(sol.x(neq+1:nw,1),[settings.nu,settings.N]);
        
    mem.mu_u_new(:) = sol.y(neq+1:neq+N*nu,1);
    mu_x = sol.y(neq+N*nu+1:neq+N*nu+N*nbx,1);
    mu = sol.y(neq+N*nu+N*nbx+1:end,1);
        
    for i=0:N-1        
        mem.mu_x_new(i*Nc*settings.nbx+1:(i+1)*Nc*settings.nbx-settings.nbx) = mu(i*nc+1:i*nc+(Nc-1)*nbx);
        mem.mu_x_new((i+1)*Nc*settings.nbx-settings.nbx+1:(i+1)*Nc*settings.nbx) = mu_x(i*settings.nbx+1:(i+1)*settings.nbx);                           
        mem.mu_new(i*Nc*settings.nc+1:(i+1)*Nc*settings.nc) = mu(i*nc+(Nc-1)*nbx+1:(i+1)*nc);
    end    
    
    cpt_qp   = sol.info.run_time*1e3;
               
    Recover(mem, settings);
end

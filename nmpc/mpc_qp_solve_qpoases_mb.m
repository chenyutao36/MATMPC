function [cpt_qp, mem] = mpc_qp_solve_qpoases_mb(sizes,mem)

    nu=sizes.nu;
    N=sizes.N;  
    nbx=sizes.nbx;
    nc=sizes.nc;
    ncN=sizes.ncN;
    
    moving_block(mem,sizes);
                   
    if mem.warm_start==0               
        [mem.warm_start,sol,fval,exitflag,iterations,multiplier,auxOutput] = qpOASES_sequence('i',mem.Hc_r,mem.gc_r,[mem.Ccx_r;mem.Ccg_r;mem.T],...
            -inf(mem.r*nu,1),inf(mem.r*nu,1),[mem.lxc;mem.lcc;mem.lb_du],[mem.uxc;mem.ucc;mem.ub_du],mem.qpoases_opt); 
    else
        if mem.hot_start==0
            [sol,fval,exitflag,iterations,multiplier,auxOutput] = qpOASES_sequence('m',mem.warm_start,mem.Hc_r,mem.gc_r,...
                [mem.Ccx_r;mem.Ccg_r;mem.T],-inf(mem.r*nu,1),inf(mem.r*nu,1),[mem.lxc;mem.lcc;mem.lb_du],[mem.uxc;mem.ucc;mem.ub_du],mem.qpoases_opt);
        end
        if mem.hot_start==1
           [sol,fval,exitflag,iterations,multiplier,auxOutput] = qpOASES_sequence('h',mem.warm_start,mem.gc,mem.lb_du,...
               mem.ub_du,[mem.lxc;mem.lcc],[mem.uxc;mem.ucc],mem.qpoases_opt);
        end
    end
    
    %The dual solution vector contains exactly one entry per lower/upper bound as well as exactly one entry per
    %lower/upper constraints bound. Positive entries correspond to active lower (constraints) bounds, negative
    %entries to active upper (constraints) bounds and a zero entry means that both corresponding (constraints)
    %bounds are inactive.
    sol = mem.T*sol;
        
    mem.du(:) = reshape(sol, [nu,N]);        
    mem.mu_x_new(:) = - multiplier(mem.r*nu+1:mem.r*nu+N*nbx);
    mem.mu_new(:)   = - multiplier(mem.r*nu+N*nbx+1:mem.r*nu+N*nbx+N*nc+ncN);
    mem.mu_u_new(:)   = - multiplier(mem.r*nu+N*nbx+N*nc+ncN+1:end);
    cpt_qp   = auxOutput.cpuTime*1e3;
                
    Recover(mem, sizes);
end

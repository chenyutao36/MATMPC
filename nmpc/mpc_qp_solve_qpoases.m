function [cpt_qp, mem] = mpc_qp_solve_qpoases(sizes,mem,opt)

    nu=sizes.nu;
    N=sizes.N;  
    nbx=sizes.nbx;
    
    if strcmp(opt.condensing,'hpipm_full')
        lb_du = flipud(mem.lb_du);
        ub_du = flipud(mem.ub_du);
    else
        lb_du = mem.lb_du;
        ub_du = mem.ub_du;
    end
           
    if mem.warm_start==0               
        [mem.warm_start,sol,fval,exitflag,iterations,multiplier,auxOutput] = qpOASES_sequence('i',mem.Hc,mem.gc,[mem.Ccx;mem.Ccg],...
            lb_du,ub_du,[mem.lxc;mem.lcc],[mem.uxc;mem.ucc],mem.qpoases_opt); 
    else
        if mem.hot_start==0
            [sol,fval,exitflag,iterations,multiplier,auxOutput] = qpOASES_sequence('m',mem.warm_start,mem.Hc,mem.gc,...
                [mem.Ccx;mem.Ccg],lb_du,ub_du,[mem.lxc;mem.lcc],[mem.uxc;mem.ucc],mem.qpoases_opt);
        end
        if mem.hot_start==1
           [sol,fval,exitflag,iterations,multiplier,auxOutput] = qpOASES_sequence('h',mem.warm_start,mem.gc,lb_du,...
               ub_du,[mem.lxc;mem.lcc],[mem.uxc;mem.ucc],mem.qpoases_opt);
        end
    end
    
    %The dual solution vector contains exactly one entry per lower/upper bound as well as exactly one entry per
    %lower/upper constraints bound. Positive entries correspond to active lower (constraints) bounds, negative
    %entries to active upper (constraints) bounds and a zero entry means that both corresponding (constraints)
    %bounds are inactive.
        
    mem.du = reshape(sol, [nu,N]);
    mem.mu_u_new  = - multiplier(1:N*nu); 
    mem.mu_x_new = -multiplier(N*nu+1:N*nu+(N+1)*nbx);
    mem.mu_new   = - multiplier(N*nu+(N+1)*nbx+1:end);
    cpt_qp   = auxOutput.cpuTime*1e3;
    
    if strcmp(opt.condensing,'hpipm_full')
        mem.du = fliplr(mem.du);
        mem.mu_u_new = flipud(mem.mu_u_new);
        mem.mu_x_new = flipud(mem.mu_x_new);
        mem.mu_new = flipud(mem.mu_new);
    end
            
    Recover(mem, sizes);
end

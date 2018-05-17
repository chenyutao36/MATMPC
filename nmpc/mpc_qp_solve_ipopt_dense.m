function [cpt_qp, mem] = mpc_qp_solve_ipopt_dense(sizes,mem)
    N=sizes.N;  
    nu=sizes.nu;
    nbx=sizes.nbx;

    H = sparse(mem.Hc);
    g = mem.gc;
   
    mem.ipopt.funcs.hessian = @(x,sigma,lambda) tril(sigma*H);
    mem.ipopt.funcs.hessianstructure = @() tril(H);
    mem.ipopt.funcs.objective=@(x) 0.5*x'*H*x + g'*x;
    mem.ipopt.funcs.gradient = @(x) (H*x + g)';
    mem.ipopt.options.A= sparse([mem.Ccx;mem.Ccg]);
    mem.ipopt.options.ru= [mem.uxc;mem.ucc];
    mem.ipopt.options.rl= [mem.lxc;mem.lcc];
    mem.ipopt.options.ub = mem.ub_du;
    mem.ipopt.options.lb = mem.lb_du;
    [mem.du,fval,exitflag,info_ipopt] = opti_ipopt(mem.ipopt,[]);
         
    mem.mu_u_new = info_ipopt.Lambda.upper-info_ipopt.Lambda.lower;
    mem.mu_x_new = info_ipopt.Lambda.ineqlin(1:N*nbx);
    mem.mu_new   = info_ipopt.Lambda.ineqlin(N*nbx+1:end);
    cpt_qp   = info_ipopt.Time*1e3;
    
    Recover(mem, sizes);

end


function [cpt_qp, mem] = mpc_qp_solve_ipopt_dense(sizes,mem)
    N=sizes.N;  
    nu=sizes.nu;
    nbx=sizes.nbx;

    H = sparse(mem.Hc);
    g = mem.gc;
    
    dBu = eye(N*nu,N*nu);

    mem.ipopt.funcs.hessian = @(x,sigma,lambda) tril(sigma*H);
    mem.ipopt.funcs.hessianstructure = @() tril(H);
    mem.ipopt.funcs.objective=@(x) 0.5*x'*H*x + g'*x;
    mem.ipopt.funcs.gradient = @(x) (H*x + g)';
    mem.ipopt.options.A= sparse([dBu;mem.Ccx;mem.Ccg]);
    mem.ipopt.options.ru= [mem.ub_du;mem.uxc;mem.ucc];
    mem.ipopt.options.rl= [mem.lb_du;mem.lxc;mem.lcc];
    [mem.du,fval,exitflag,info_ipopt] = opti_ipopt(mem.ipopt,[]);
         
    mem.mu_u_new = info_ipopt.Lambda.ineqlin(1:N*nu);
    mem.mu_x_new = info_ipopt.Lambda.ineqlin(N*nu+1:N*nu+N*nbx);
    mem.mu_new   = info_ipopt.Lambda.ineqlin(N*nu+N*nbx+1:end);
    cpt_qp   = info_ipopt.Time*1e3;
    
    Recover(mem, sizes);

end


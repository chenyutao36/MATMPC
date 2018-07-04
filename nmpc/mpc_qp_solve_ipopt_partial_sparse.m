function [cpt_qp, mem] = mpc_qp_solve_ipopt_partial_sparse(settings,settings2,mem, mem2)
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
    
    ipopt.options = mem2.ipopt.options;
    ipopt.x0 = mem2.ipopt.x0;
            
    for i=0:N-1       
       for j=1:nbx
           ipopt.options.ub((i+1)*nx+nbx_idx(j))= mem2.ub_dx(i*nbx+j);
           ipopt.options.lb((i+1)*nx+nbx_idx(j))= mem2.lb_dx(i*nbx+j);
       end
    end
    
    ipopt.options.ub(neq+1:end) = mem2.ub_du;
    ipopt.options.lb(neq+1:end) = mem2.lb_du;
    
    full2sparse(mem2,settings2);
     
    H = sparse(mem2.sparse_H);
        
    ipopt.funcs.hessian = @(x,sigma,lambda) tril(sigma*H);
    ipopt.funcs.hessianstructure = @() tril(H);
    ipopt.funcs.objective=@(x) 0.5*x'*H*x + (mem2.sparse_g)'*x;
    ipopt.funcs.gradient = @(x) (H*x + mem2.sparse_g)';
    ipopt.options.A= sparse([mem2.sparse_dB;mem2.sparse_dG]);
    ipopt.options.ru= [mem2.uc;-mem2.sparse_G];
    ipopt.options.rl= [mem2.lc;-mem2.sparse_G];
    [sol,fval,exitflag,info_ipopt] = opti_ipopt(ipopt,[]);
    cpt_qp   = info_ipopt.Time*1e3;
    
    mem.du(:) = reshape(sol(neq+1:nw,1),[settings.nu,settings.N]);
    mem.mu_u_new(:)   = info_ipopt.Lambda.upper((N+1)*nx+1:end) - info_ipopt.Lambda.lower((N+1)*nx+1:end);
    mu_x_p = info_ipopt.Lambda.upper(1:(N+1)*nx) - info_ipopt.Lambda.lower(1:(N+1)*nx);
    mu     = info_ipopt.Lambda.ineqlin;
    
    mu_x   = zeros(N*nbx,1);
    for i=0:N-1
        for j=1:nbx
            mu_x(i*nbx+j,1)  = mu_x_p((i+1)*nx+nbx_idx(j));
        end
    end
    
    for i=0:N-1        
        mem.mu_x_new(i*Nc*settings.nbx+1:(i+1)*Nc*settings.nbx-settings.nbx) = mu(i*nc+1:i*nc+(Nc-1)*nbx);
        mem.mu_x_new((i+1)*Nc*settings.nbx-settings.nbx+1:(i+1)*Nc*settings.nbx) = mu_x(i*settings.nbx+1:(i+1)*settings.nbx);                   
    end
    
    Recover(mem, settings);
end


function [cpt_qp, mem] = mpc_qp_solve_ipopt_sparse(sizes,mem)
    N=sizes.N;  
    nbx=sizes.nbx;
    nbx_idx=sizes.nbx_idx;
    nx=sizes.nx;
    nu=sizes.nu;
    nc=sizes.nc;
    ncN=sizes.ncN;
    neq = (N+1)*nx;
    nw = (N+1)*nx+N*nu;
    
    ipopt.options = mem.ipopt.options;
    ipopt.x0 = mem.ipopt.x0;

    for i=0:N-1      
       for j=1:nbx
           ipopt.options.ub((i+1)*nx+nbx_idx(j))= mem.ub_dx(i*nbx+j);
           ipopt.options.lb((i+1)*nx+nbx_idx(j))= mem.lb_dx(i*nbx+j);
       end
    end

    ipopt.options.ub(neq+1:end) = mem.ub_du;
    ipopt.options.lb(neq+1:end) = mem.lb_du;

    full2sparse(mem,sizes);
    
    H = sparse(mem.sparse_H);
    
    ipopt.funcs.hessian = @(x,sigma,lambda) tril(sigma*H);
    ipopt.funcs.hessianstructure = @() tril(H);
    ipopt.funcs.objective=@(x) 0.5*x'*H*x + (mem.sparse_g)'*x;
    ipopt.funcs.gradient = @(x) (H*x + mem.sparse_g)';
    ipopt.options.A= sparse([mem.sparse_dB;mem.sparse_dG]);
    ipopt.options.ru= [mem.uc;-mem.sparse_G];
    ipopt.options.rl= [mem.lc;-mem.sparse_G];
    [sol,fval,exitflag,info_ipopt] = opti_ipopt(ipopt,[]);
            
    mem.dx=reshape(sol(1:neq,1),[nx,N+1]);
    mem.du=reshape(sol(neq+1:nw,1),[nu,N]);
            
    mem.mu_u_new   = info_ipopt.Lambda.upper((N+1)*nx+1:end) - info_ipopt.Lambda.lower((N+1)*nx+1:end);
    mu_x = info_ipopt.Lambda.upper(1:(N+1)*nx) - info_ipopt.Lambda.lower(1:(N+1)*nx);
    for i=0:N-1
        for j=1:nbx
            mem.mu_x_new(i*nbx+j,1)  = mu_x((i+1)*nx+nbx_idx(j));
        end
    end
    mem.mu_new     = info_ipopt.Lambda.ineqlin;
    mem.lambda_new = reshape(info_ipopt.Lambda.eqlin,[nx N+1]);
    cpt_qp   = info_ipopt.Time*1e3;
end


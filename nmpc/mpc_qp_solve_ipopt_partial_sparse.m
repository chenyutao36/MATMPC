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
            
    mem2.ipopt_data.G(1:nx,1) = mem2.ds0(:);
    for i=0:N-1
       mem2.ipopt_data.H(i*nx+1:(i+1)*nx, i*nx+1:(i+1)*nx)  = mem2.Q(:,i*nx+1:(i+1)*nx);
       mem2.ipopt_data.H(i*nx+1:(i+1)*nx, neq+i*nu+1:neq+(i+1)*nu)  = mem2.S(:,i*nu+1:(i+1)*nu);
       mem2.ipopt_data.H(neq+i*nu+1:neq+(i+1)*nu, i*nx+1:(i+1)*nx)  = (mem2.S(:,i*nu+1:(i+1)*nu))';
       mem2.ipopt_data.H(neq+i*nu+1:neq+(i+1)*nu, neq+i*nu+1:neq+(i+1)*nu)  = mem2.R(:,i*nu+1:(i+1)*nu);
       
       mem2.ipopt_data.dBg(i*nc+1:(i+1)*nc, i*nx+1:(i+1)*nx) = mem2.Cgx(:,i*nx+1:(i+1)*nx);
       mem2.ipopt_data.dBg(i*nc+1:(i+1)*nc, neq+i*nu+1:neq+(i+1)*nu) = mem2.Cgu(:,i*nu+1:(i+1)*nu);
       
       mem2.ipopt_data.dG((i+1)*nx+1:(i+2)*nx, i*nx+1:(i+2)*nx) = [mem2.A(:,i*nx+1:(i+1)*nx), -eye(nx,nx)];
       mem2.ipopt_data.dG((i+1)*nx+1:(i+2)*nx, neq+i*nu+1:neq+(i+1)*nu) = mem2.B(:,i*nu+1:(i+1)*nu);
       
       mem2.ipopt_data.G((i+1)*nx+1:(i+2)*nx, 1) = -mem2.a(:,i+1);
       
       mem2.ipopt_data.g(i*nx+1:(i+1)*nx) = mem2.gx(:,i+1);
       mem2.ipopt_data.g(neq+i*nu+1:neq+(i+1)*nu) = mem2.gu(:,i+1);
       
       for j=1:nbx
           ipopt.options.ub((i+1)*nx+nbx_idx(j))= mem2.ub_dx(i*nbx+j);
           ipopt.options.lb((i+1)*nx+nbx_idx(j))= mem2.lb_dx(i*nbx+j);
       end
    end
    mem2.ipopt_data.H(N*nx+1:(N+1)*nx, N*nx+1:(N+1)*nx)  = mem2.Q(:,N*nx+1:(N+1)*nx);
    mem2.ipopt_data.g(N*nx+1:(N+1)*nx) = mem2.gx(:,N+1);
    mem2.ipopt_data.dBg(N*nc+1:N*nc+ncN, N*nx+1:(N+1)*nx) = mem2.CgN;
    
    ipopt.options.ub(neq+1:end) = mem2.ub_du;
    ipopt.options.lb(neq+1:end) = mem2.lb_du;
     
    H = sparse(mem2.ipopt_data.H);
        
    ipopt.funcs.hessian = @(x,sigma,lambda) tril(sigma*H);
    ipopt.funcs.hessianstructure = @() tril(H);
    ipopt.funcs.objective=@(x) 0.5*x'*H*x + (mem2.ipopt_data.g)'*x;
    ipopt.funcs.gradient = @(x) (H*x + mem2.ipopt_data.g)';
    ipopt.options.A= sparse([mem2.ipopt_data.dBg;mem2.ipopt_data.dG]);
    ipopt.options.ru= [mem2.uc;mem2.ipopt_data.G];
    ipopt.options.rl= [mem2.lc;mem2.ipopt_data.G];
    [sol,fval,exitflag,info_ipopt] = opti_ipopt(ipopt,[]);
    cpt_qp   = info_ipopt.Time*1e3;
    
    mem.du = reshape(sol(neq+1:nw,1),[settings.nu,settings.N]);
    mem.mu_u_new   = info_ipopt.Lambda.upper((N+1)*nx+1:end) - info_ipopt.Lambda.lower((N+1)*nx+1:end);
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


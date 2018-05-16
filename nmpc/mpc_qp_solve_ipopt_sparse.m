function [cpt_qp, mem] = mpc_qp_solve_ipopt_sparse(sizes,mem)
    N=sizes.N;  
    nbx=sizes.nbx;
    nx=sizes.nx;
    nu=sizes.nu;
    nc=sizes.nc;
    ncN=sizes.ncN;
    neq = (N+1)*nx;
    nw = (N+1)*nx+N*nu;
    
    mem.ipopt_data.G(1:nx,1) = mem.ds0(:);
    for i=0:N-1
       mem.ipopt_data.H(i*nx+1:(i+1)*nx, i*nx+1:(i+1)*nx)  = mem.Q(:,i*nx+1:(i+1)*nx);
       mem.ipopt_data.H(neq+i*nu+1:neq+(i+1)*nu, neq+i*nu+1:neq+(i+1)*nu)  = mem.R(:,i*nu+1:(i+1)*nu);
       
       mem.ipopt_data.dBx(i*nbx+1:(i+1)*nbx, (i+1)*nx+1:(i+2)*nx) = mem.Cx;
       mem.ipopt_data.dBg(i*nc+1:(i+1)*nc, i*nx+1:(i+1)*nx) = mem.Cgx(:,i*nx+1:(i+1)*nx);
       mem.ipopt_data.dBg(i*nc+1:(i+1)*nc, neq+i*nu+1:neq+(i+1)*nu) = mem.Cgu(:,i*nu+1:(i+1)*nu);
       
       mem.ipopt_data.dG((i+1)*nx+1:(i+2)*nx, i*nx+1:(i+2)*nx) = [mem.A(:,i*nx+1:(i+1)*nx), -eye(nx,nx)];
       mem.ipopt_data.dG((i+1)*nx+1:(i+2)*nx, neq+i*nu+1:neq+(i+1)*nu) = mem.B(:,i*nu+1:(i+1)*nu);
       
       mem.ipopt_data.G((i+1)*nx+1:(i+2)*nx, 1) = -mem.a(:,i+1);
       
       mem.ipopt_data.g(i*nx+1:(i+1)*nx) = mem.gx(:,i+1);
       mem.ipopt_data.g(neq+i*nu+1:neq+(i+1)*nu) = mem.gu(:,i+1);
    end
    mem.ipopt_data.H(N*nx+1:(N+1)*nx, N*nx+1:(N+1)*nx)  = mem.Q(:,N*nx+1:(N+1)*nx);
    mem.ipopt_data.g(N*nx+1:(N+1)*nx) = mem.gx(:,N+1);
    mem.ipopt_data.dBg(N*nc+1:N*nc+ncN, N*nx+1:(N+1)*nx) = mem.CgN;
    
    H = sparse(mem.ipopt_data.H);
    
    mem.ipopt.funcs.hessian = @(x,sigma,lambda) tril(sigma*H);
    mem.ipopt.funcs.hessianstructure = @() tril(H);
    mem.ipopt.funcs.objective=@(x) 0.5*x'*H*x + (mem.ipopt_data.g)'*x;
    mem.ipopt.funcs.gradient = @(x) (H*x + mem.ipopt_data.g)';
    mem.ipopt.options.A= sparse([mem.ipopt_data.dBu;mem.ipopt_data.dBx;mem.ipopt_data.dBg;mem.ipopt_data.dG]);
    mem.ipopt.options.ru= [mem.ub_du;mem.ub_dx;mem.uc;mem.ipopt_data.G];
    mem.ipopt.options.rl= [mem.lb_du;mem.lb_dx;mem.lc;mem.ipopt_data.G];
    [sol,fval,exitflag,info_ipopt] = opti_ipopt(mem.ipopt,[]);
            
    mem.dx=reshape(sol(1:neq,1),[nx,N+1]);
    mem.du=reshape(sol(neq+1:nw,1),[nu,N]);
            
    mem.mu_u_new   = info_ipopt.Lambda.ineqlin(1:N*nu);
    mem.mu_x_new   = info_ipopt.Lambda.ineqlin(N*nu+1:N*nu+N*nbx);
    mem.mu_new     = info_ipopt.Lambda.ineqlin(N*nu+N*nbx+1:end);
    mem.lambda_new = reshape(info_ipopt.Lambda.eqlin,[nx N+1]);
    cpt_qp   = info_ipopt.Time*1e3;
end


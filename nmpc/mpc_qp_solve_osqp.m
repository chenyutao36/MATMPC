function [cpt_qp, mem] = mpc_qp_solve_osqp(sizes,mem)

    N=sizes.N;  
    nbx=sizes.nbx;
    nx=sizes.nx;
    nu=sizes.nu;
    nc=sizes.nc;
    ncN=sizes.ncN;
    neq = (N+1)*nx;
    nw = (N+1)*nx+N*nu;
    
    mem.osqp_data.G(1:nx,1) = mem.ds0(:);
    
    for i=0:N-1
       mem.osqp_data.H(i*nx+1:(i+1)*nx, i*nx+1:(i+1)*nx)  = mem.Q(:,i*nx+1:(i+1)*nx);
       mem.osqp_data.H(i*nx+1:(i+1)*nx, neq+i*nu+1:neq+(i+1)*nu)  = mem.S(:,i*nu+1:(i+1)*nu);
       mem.osqp_data.H(neq+i*nu+1:neq+(i+1)*nu, i*nx+1:(i+1)*nx)  = (mem.S(:,i*nu+1:(i+1)*nu))';
       mem.osqp_data.H(neq+i*nu+1:neq+(i+1)*nu, neq+i*nu+1:neq+(i+1)*nu)  = mem.R(:,i*nu+1:(i+1)*nu);
       
       mem.osqp_data.dBg(i*nc+1:(i+1)*nc, i*nx+1:(i+1)*nx) = mem.Cgx(:,i*nx+1:(i+1)*nx);
       mem.osqp_data.dBg(i*nc+1:(i+1)*nc, neq+i*nu+1:neq+(i+1)*nu) = mem.Cgu(:,i*nu+1:(i+1)*nu);
       
       mem.osqp_data.dG((i+1)*nx+1:(i+2)*nx, i*nx+1:(i+2)*nx) = [mem.A(:,i*nx+1:(i+1)*nx), -eye(nx,nx)];
       mem.osqp_data.dG((i+1)*nx+1:(i+2)*nx, neq+i*nu+1:neq+(i+1)*nu) = mem.B(:,i*nu+1:(i+1)*nu);
       
       mem.osqp_data.G((i+1)*nx+1:(i+2)*nx, 1) = -mem.a(:,i+1);
       
       mem.osqp_data.g(i*nx+1:(i+1)*nx) = mem.gx(:,i+1);
       mem.osqp_data.g(neq+i*nu+1:neq+(i+1)*nu) = mem.gu(:,i+1);
       
    end
    mem.osqp_data.H(N*nx+1:(N+1)*nx, N*nx+1:(N+1)*nx)  = mem.Q(:,N*nx+1:(N+1)*nx);
    mem.osqp_data.g(N*nx+1:(N+1)*nx) = mem.gx(:,N+1);
    mem.osqp_data.dBg(N*nc+1:N*nc+ncN, N*nx+1:(N+1)*nx) = mem.CgN;
    
    mem.osqp_data.ub(1:N*nu) = mem.ub_du;
    mem.osqp_data.lb(1:N*nu) = mem.lb_du;
    mem.osqp_data.ub(N*nu+1:N*nu+N*nbx) = mem.ub_dx;
    mem.osqp_data.lb(N*nu+1:N*nu+N*nbx) = mem.lb_dx;
    mem.osqp_data.ub(N*nu+N*nbx+1:end) = mem.uc;
    mem.osqp_data.lb(N*nu+N*nbx+1:end) = mem.lc;
    
    A = [mem.osqp_data.dG; mem.osqp_data.dBu; mem.osqp_data.dBx; mem.osqp_data.dBg];
    l = [mem.osqp_data.G; mem.osqp_data.lb];
    u = [mem.osqp_data.G; mem.osqp_data.ub];
    
%     if mem.iter==1 && mem.sqp_it== 0               
%         mem.qp_obj.setup(mem.osqp_data.H, mem.osqp_data.g, A, l, u, mem.osqp_options);
%         sol = mem.qp_obj.solve();
%     else
%         mem.qp_obj.update('Px', mem.osqp_data.H, 'Ax', A);
%         mem.qp_obj.update('q', mem.osqp_data.g, 'l', l, 'u', u);
%         sol = mem.qp_obj.solve();
%     end

    qp_obj = osqp;
    qp_obj.setup(mem.osqp_data.H, mem.osqp_data.g, A, l, u, mem.osqp_options);
    sol = qp_obj.solve();       

    mem.dx(:) = reshape(sol.x(1:neq,1),[nx,N+1]);
    mem.du(:) = reshape(sol.x(neq+1:nw,1),[nu,N]);
        
    mem.lambda_new(:) = reshape(sol.y(1:neq,1),[nx,N+1]);
    mem.mu_u_new(:) = sol.y(neq+1:neq+N*nu,1);
    mem.mu_x_new(:) = sol.y(neq+N*nu+1:neq+N*nu+N*nbx,1);
    mem.mu_new(:) = sol.y(neq+N*nu+N*nbx+1:end,1);
    
    cpt_qp   = sol.info.run_time*1e3;
               
end

function [cpt_qp, mem] = mpc_qp_solve_qpdunes(sizes,mem)

    nx=sizes.nx;
    nu=sizes.nu;
    N=sizes.N;  
    nbx=sizes.nbx;
    nbx_idx=sizes.nbx_idx;
    nc=sizes.nc;
    ncN=sizes.ncN;
    
    for i=0:N-1
        mem.dunes.H(1:nx,i*(nx+nu)+1:i*(nx+nu)+nx) = mem.Q(:,i*nx+1:(i+1)*nx);
        mem.dunes.H(1:nx,i*(nx+nu)+nx+1:(i+1)*(nx+nu)) = mem.S(:,i*nu+1:(i+1)*nu);
        mem.dunes.H(nx+1:nx+nu,i*(nx+nu)+1:i*(nx+nu)+nx) = (mem.S(:,i*nu+1:(i+1)*nu))';
        mem.dunes.H(nx+1:nx+nu,i*(nx+nu)+nx+1:(i+1)*(nx+nu)) = mem.R(:,i*nu+1:(i+1)*nu);
        
        mem.dunes.C(:,i*(nx+nu)+1:(i+1)*(nx+nu)) = [mem.A(:,i*nx+1:(i+1)*nx),mem.B(:,i*nu+1:(i+1)*nu)];
        mem.dunes.c(:,i+1) = mem.a(:,i+1);   
        
        mem.dunes.D(:,i*(nx+nu)+1:(i+1)*(nx+nu)) =[mem.Cgx(:,i*nx+1:(i+1)*nx),mem.Cgu(:,i*nu+1:(i+1)*nu)];
        
        mem.dunes.g(i*(nx+nu)+1:(i+1)*(nx+nu)) = [mem.gx(:,i+1);mem.gu(:,i+1)];
    end
    mem.dunes.P = mem.Q(:,N*nx+1:(N+1)*nx);
    mem.dunes.g(N*(nx+nu)+1:end) = mem.gx(:,N+1);
    mem.dunes.D(1:ncN,N*(nx+nu)+1:N*(nx+nu)+nx) = mem.CgN;
        
    mem.dunes.zUpp(1:nx,1)=mem.ds0;
    mem.dunes.zLow(1:nx,1)=mem.ds0;
    mem.dunes.zUpp(nx+1:nx+nu,1)=mem.ub_du(1:nu,1);
    mem.dunes.zLow(nx+1:nx+nu,1)=mem.lb_du(1:nu,1);
    for i=1:N-1
        for j=1:nbx
           mem.dunes.zUpp(i*(nx+nu)+nbx_idx(j))= mem.ub_dx((i-1)*nbx+j);
           mem.dunes.zLow(i*(nx+nu)+nbx_idx(j))= mem.lb_dx((i-1)*nbx+j);
        end
        mem.dunes.zUpp(i*(nx+nu)+nx+1:(i+1)*(nx+nu),1)=mem.ub_du(i*nu+1:(i+1)*nu,1);
        mem.dunes.zLow(i*(nx+nu)+nx+1:(i+1)*(nx+nu),1)=mem.lb_du(i*nu+1:(i+1)*nu,1);
    end
    for j=1:nbx
       mem.dunes.zUpp(N*(nx+nu)+nbx_idx(j))= mem.ub_dx((N-1)*nbx+j);
       mem.dunes.zLow(N*(nx+nu)+nbx_idx(j))= mem.lb_dx((N-1)*nbx+j);
    end
    
    mem.dunes.dUpp(1:N*nc+ncN) = mem.uc;
    mem.dunes.dLow(1:N*nc+ncN) = mem.lc;
    
    tic;
               
    if mem.iter==1               
%         qpDUNES( 'init', N, ...
%           mem.dunes.H, mem.dunes.P, mem.dunes.g, ...
%           mem.dunes.C, mem.dunes.c, ...
%           mem.dunes.zLow, mem.dunes.zUpp, ...
%           [], [], [], ...       
%           mem.dunes.qpOptions );
      
      qpDUNES_aff( 'init', N, ...
          mem.dunes.H, mem.dunes.P, mem.dunes.g, ...
          mem.dunes.C, mem.dunes.c, ...
          mem.dunes.zLow, mem.dunes.zUpp, ...
          mem.dunes.D, mem.dunes.dLow, mem.dunes.dUpp, ...       
          mem.dunes.qpOptions );
    else
%         qpDUNES( 'update', mem.dunes.H, mem.dunes.P, mem.dunes.g, mem.dunes.C, mem.dunes.c, mem.dunes.zLow, mem.dunes.zUpp, [], [], [] );
        qpDUNES_aff( 'update', mem.dunes.H, mem.dunes.P, mem.dunes.g, mem.dunes.C, mem.dunes.c, mem.dunes.zLow, mem.dunes.zUpp, mem.dunes.D, mem.dunes.dLow, mem.dunes.dUpp);
    end
    
%     [sol, stat, lambda, multiplier, objFctnVal] = qpDUNES( 'solve' );
    [sol, stat, lambda, multiplier, objFctnVal] = qpDUNES_aff( 'solve' );
    
    cpt_qp   = toc*1e3;
     
    for i=0:N-1
        mem.du(:,i+1) = sol(i*(nx+nu)+nx+1:(i+1)*(nx+nu));
        mem.dx(:,i+1) = sol(i*(nx+nu)+1:i*(nx+nu)+nx);
    end
    mem.dx(:,N+1) = sol(N*(nx+nu)+1:N*(nx+nu)+nx);
    
    mem.lambda_new = reshape(lambda,[nx,N]);
        
    if ~isempty(multiplier)
        mem.mu_u_new  = - multiplier(1:N*nu); 
        mem.mu_x_new = -multiplier(N*nu+1:N*nu+N*nbx);
        mem.mu_new   = - multiplier(N*nu+N*nbx+1:end);
    end
    
              
end

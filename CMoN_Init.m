function [rho, c2] = CMoN_Init(sizes, mem, input)

    nx=sizes.nx;
    nu=sizes.nu;
    nc=sizes.nc;
    ncN=sizes.ncN;
    N=sizes.N;
    nbu=sizes.nbu;
    nbu_idx=sizes.nbu_idx;
    nbx=sizes.nbx;
    nbx_idx=sizes.nbx_idx;

    nz = nx+nu;
    nw = (N+1)*nx+N*nu;
    neq = (N+1)*nx;
    nineq = 2*(N*nc+ncN+N*nbu+N*nbx);

    H=[];
    for i=0:N-1
        H=blkdiag(H, [mem.Q(:,i*nx+1:(i+1)*nx), mem.S(:,i*nu+1:(i+1)*nu); 
                      (mem.S(:,i*nu+1:(i+1)*nu))', mem.R(:,i*nu+1:(i+1)*nu)]);    
    end
    H=blkdiag(H, mem.Q(:,N*nx+1:(N+1)*nx));

    BT=zeros(nw,neq);
    BT(1:nz,1:2*nx)=[eye(nx,nx), (mem.A(:,1:nx))'; zeros(nu,nx), (mem.B(:,1:nu))'];
    for i=1:N-1
        BT(i*nz+1:(i+1)*nz, i*nx+1:(i+2)*nx) = [-eye(nx,nx), (mem.A(:,i*nx+1:(i+1)*nx))'; zeros(nu,nx), (mem.B(:,i*nu+1:(i+1)*nu))'];
    end
    BT(N*nz+1:N*nz+nx, N*nx+1:(N+1)*nx) = -eye(nx,nx);  
    B=BT';

    C=[];
    Iu = zeros(nbu,nu);
    for j=1:nbu
        Iu(j,nbu_idx(j))=1;
    end
    Ix = zeros(nbx,nx);
    for j=1:nbx
        Ix(j,nbx_idx(j))=1;
    end
    for i=0:N-1
        if i>0
            C=blkdiag(C, [mem.Cgx(:,i*nx+1:(i+1)*nx), mem.Cgu(:,i*nu+1:(i+1)*nu);
                          -mem.Cgx(:,i*nx+1:(i+1)*nx) , -mem.Cgu(:,i*nu+1:(i+1)*nu);
                          zeros(nbu,nx), Iu;
                          zeros(nbu,nx), -Iu;
                          Ix, zeros(nbx,nu);
                          -Ix, zeros(nbx,nu)]);
        else
            C=blkdiag(C, [mem.Cgx(:,i*nx+1:(i+1)*nx), mem.Cgu(:,i*nu+1:(i+1)*nu);
                          -mem.Cgx(:,i*nx+1:(i+1)*nx) , -mem.Cgu(:,i*nu+1:(i+1)*nu);
                          zeros(nbu,nx), Iu;
                          zeros(nbu,nx), -Iu]);
        end
    end
    C=blkdiag(C,[mem.CgN;-mem.CgN;
                 Ix;
                 -Ix]);

    blk4 = [];    
    no_mu = zeros(1,N+1);
    no_mu(1) = 2*nc+2*nbu;
    for i=2:N
        no_mu(i) = 2*(nc+nbu+nbx);
    end
    no_mu(N+1) = 2*(ncN+nbx); 
    
    dmu = mem.dmu(N*nu+N*nbx+1:end);
    dmu_u = mem.dmu(1:N*nu);
    dmu_x= mem.dmu(N*nu+1:N*nu+N*nbx);
%     dmu = mem.mu_new;
%     dmu_u = mem.mu_u_new;
%     dmu_x = mem.mu_x_new;

    for i=0:N-1
        if i>0
            mu_tmp = [max(dmu(i*nc+1:(i+1)*nc),0);min(dmu(i*nc+1:(i+1)*nc),0);
                      max(dmu_u(i*nbu+1:(i+1)*nbu),0);min(dmu_u(i*nbu+1:(i+1)*nbu),0);
                      max(dmu_x((i-1)*nbx+1:(i)*nbx),0);min(dmu_x((i-1)*nbx+1:(i)*nbx),0)];
            blk4=[blk4; -mu_tmp.*C(sum(no_mu(1:i))+1:sum(no_mu(1:i+1)),:)];
        else
            mu_tmp = [max(dmu(i*nc+1:(i+1)*nc),0);min(dmu(i*nc+1:(i+1)*nc),0);
                      max(dmu_u(i*nbu+1:(i+1)*nbu),0);min(dmu_u(i*nbu+1:(i+1)*nbu),0)];
            blk4=[blk4; -mu_tmp.*C(1:no_mu(1),:)];
        end
      
    end
    mu_tmp = [max(dmu(N*nc+1:N*nc+ncN),0);min(dmu(N*nc+1:N*nc+ncN),0);
              max(dmu_x((N-1)*nbx+1:N*nbx),0);min(dmu_x((N-1)*nbx+1:N*nbx),0)];
    blk4=[blk4; -mu_tmp.*C(sum(no_mu(1:N))+1:sum(no_mu(1:N+1)),:)];
  
    c=[];
    tmp1=[];
    tmp2=[];
    for i=0:N-1
        if i>0
            c=[c;-mem.uc(i*nc+1:(i+1)*nc);
                 mem.lc(i*nc+1:(i+1)*nc)];                 
                 if nbu_idx>0
                    tmp1 = mem.ub_du(i*nu+1:(i+1)*nu);
                    tmp2 = mem.lb_du(i*nu+1:(i+1)*nu);
                    tmp1 = tmp1(nbu_idx);
                    tmp2 = tmp2(nbu_idx);
                 end            
                c=[c;-tmp1;
                     tmp2;
                     -mem.ub_dx((i-1)*nbx+1:i*nbx);
                     mem.lb_dx((i-1)*nbx+1:i*nbx)];
        else
            c=[c;-mem.uc(i*nc+1:(i+1)*nc);
                 mem.lc(i*nc+1:(i+1)*nc)];                 
                 if nbu_idx>0
                    tmp1 = mem.ub_du(i*nu+1:(i+1)*nu);
                    tmp2 = mem.lb_du(i*nu+1:(i+1)*nu);
                    tmp1 = tmp1(nbu_idx);
                    tmp2 = tmp2(nbu_idx);
                 end        
            c=[c;-tmp1;
                 tmp2;];
        end
    end
    c=[c;-mem.uc(N*nc+1:end);
         mem.lc(N*nc+1:end);
         -mem.ub_dx((N-1)*nbx+1:N*nbx);
         mem.lb_dx((N-1)*nbx+1:N*nbx)];
    blk5 = -diag(c);
   
    M=[H,C',BT;
       blk4, blk5, zeros(nineq, neq);
       B, zeros(neq,nineq), zeros(neq, neq)];
   
    s = svd(M).^-1;
    rho = max(s);   
    c2 = std(s) +1;

end

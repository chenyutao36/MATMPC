function [dw,lambda,mu,info,opt] = mpc_qp_condense_and_solve(sizes,x0,H_mat,g_vec,G_vec,dG_mat,B_vec,dB_mat,z,input)

    N=sizes.N;
    nx=sizes.nx;
    nu=sizes.nu;
    nz=nx+nu;
    
    switch (input.opt.condensing)
        case 'null-space'            
                    
            tcond=tic; % time for null-space factorization
            
            [Y,R,e]=qr(dG_mat',0); 
            L = R';
            Z = spnull(Y');
            Hc   = Z'*H_mat*Z;
            sy   = -mldivide(L,G_vec(e,1));
            D    = Y'*H_mat;
            gc   = Z'*(g_vec+D'*sy);
            Bc   = B_vec+dB_mat*Y*sy;
            dBc  = dB_mat*Z;
            
            tCOND=toc(tcond)*1e3;
            
            [sz,~,mu,~,info,opt] = mpc_qp_solve(sizes,Hc,gc,dBc,Bc,input.opt,[]);

            dw=Y*sy+Z*sz;
            
            lambda_pm=mldivide(R, -Y'*dB_mat'*mu-Y'*g_vec-D*dw);
            lambda(e,1)=lambda_pm;  % reverse the permutation

        case 'standard'
                       
            mG=input.opt.condensing_matrix.mG;
            mG(1:nx,1)= -z(1,1:nx)';
            
            M=input.opt.condensing_matrix.M;
            L=input.opt.condensing_matrix.L;
            
            dG_full=full(dG_mat);           
            
            tcond=tic;
            for i=1:N

                Gxk=dG_full(i*nx+1:(i+1)*nx,(i-1)*nz+1:i*nz-nu);
                Guk=dG_full(i*nx+1:(i+1)*nx,(i-1)*nz+nx+1:i*nz);

                M(i*nz+1:(i+1)*nz-nu,1:i*nu)=Gxk*M((i-1)*nz+1:i*nz-nu,1:i*nu);
                M(i*nz+1:(i+1)*nz-nu,(i-1)*nu+1:i*nu)=Guk;
                
                if i<N
                    M(i*nz+nx+1:(i+1)*nz,i*nu+1:(i+1)*nu)=eye(nu);
                end
 
                L(i*nz+1:(i+1)*nz-nu,:) = Gxk*L((i-1)*nz+1:i*nz-nu,:);      
                mG(i*nz+1:(i+1)*nz-nu,:) = Gxk*mG((i-1)*nz+1:i*nz-nu,:)+G_vec(i*nx+1:(i+1)*nx,1);
                                             
            end
            
            tmp1 = M'*H_mat;           
        	Hc = tmp1*M;
            tmp2 = L*x0;
            tmp3 = mG+tmp2;
            gc = tmp1*tmp3+M'*g_vec;
            Bc = B_vec+dB_mat*tmp3;
        	dBc = dB_mat*M;

            tCOND=toc(tcond)*1e3;
            
            [du,~,mu,~,info,opt] = mpc_qp_solve(sizes,Hc,gc,dBc,Bc,input.opt,[]);

            dw = mG+ tmp2+ M*du;
            
            % recover the multiplier of the equality constraints
            bKKT = -H_mat*dw-g_vec-dB_mat'*mu;
            lambda = mldivide(dG_mat(:,input.data.index_lam)',bKKT(input.data.index_lam));
            
        %% end case 'standard'
    end %% end switch
    
    info.tcond=tCOND;
end
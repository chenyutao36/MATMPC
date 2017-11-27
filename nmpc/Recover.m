function [ dz,dxN, lambda,mu,muN ] = Recover( Q,S,A,B,Cx,a,gx,du,ds0,mu_vec,input,sizes )

% N=sizes.N;
% nx=sizes.nx;
% nu=sizes.nu;
% nc=sizes.nc*2;
% ncN=sizes.ncN*2;

% mu=reshape(mu_vec(1:N*nc), [nc N]);
% muN=mu_vec(N*nc+1:N*nc+ncN);

% lambda=input.lambda;
% dw=zeros((N+1)*nx+nu,1);
% 
% dw(1:nx,1)=ds0;
% for i=1:N
%      dw(i*nz+1:(i+1)*nz-nu)=A{i}*dw((i-1)*nz+1:i*nz-nu)+B{i}*du((i-1)*nu+1:i*nu,1)+a(:,i);
%      dw((i-1)*nz+nx+1:i*nz)=du((i-1)*nu+1:i*nu,1);
% end
%         
% lambda(:,N+1) = Q{N+1}*dw(N*nz+1:(N+1)*nz-nu)+gx(:,N+1)+Cx{N+1}'*muN;
% for i=N:-1:1
%      lambda(:,i)= A{i}'*lambda(:,i+1)+Q{i}*dw((i-1)*nz+1:i*nz-nu)+gx(:,i)+S{i}*du((i-1)*nu+1:i*nu,1)+Cx{i}'*mu(:,i);
% end

% dz=reshape(dw(1:N*(nx+nu)), [nx+nu N]);
% dxN=dw(N*(nx+nu)+1:end,1);  

[dz, dxN, lambda, mu, muN] = Recover_mex(Q,S,A,B,Cx,a,gx,du,ds0,mu_vec,sizes);

end


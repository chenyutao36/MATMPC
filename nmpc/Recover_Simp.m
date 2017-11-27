function [dw] = Recover_Simp( A,B,a,du,ds0,sizes )

N=sizes.N;
nx=sizes.nx;
nu=sizes.nu;

nz=nx+nu;

dw=zeros((N+1)*nx+nu,1);

dw(1:nx,1)=ds0;
for i=1:N
    dw(i*nz+1:(i+1)*nz-nu)=A{i}*dw((i-1)*nz+1:i*nz-nu)+B{i}*du((i-1)*nu+1:i*nu,1)+a(:,i);
    dw((i-1)*nz+nx+1:i*nz)=du((i-1)*nu+1:i*nu,1);
end
        

end


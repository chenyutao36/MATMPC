function [mem2] = Pcond(mem, settings, mem2, settings2)
    
    nbx = settings.nbx;
    N = settings.N;

    nx = settings2.nx;
    nu = settings2.nu;
    N2 = settings2.N;
    
    Npc = N/N2;
    
    [Hp,gp,Ccx,Ap,Bp,ap,lxc,uxc]=partial_condensing_default(mem, settings);
    
    mem2.A=Ap;
    mem2.B=Bp;
    mem2.a=ap;
    mem2.lb_du=mem.lb_du;
    mem2.ub_du=mem.ub_du;    
    mem2.lc=lxc;
    mem2.uc=uxc;
    mem2.ds0=mem.ds0;
    mem2.CgN=mem.CgN;
                                
    for i=0:N2-1
        mem2.Q(:,i*nx+1:(i+1)*nx) = Hp(1:nx,i*(nx+nu)+1:i*(nx+nu)+nx);       
        mem2.S(:,i*nu+1:(i+1)*nu) = Hp(1:nx,i*(nx+nu)+nx+1:(i+1)*(nx+nu));        
        mem2.R(:,i*nu+1:(i+1)*nu) = Hp(nx+1:end,i*(nx+nu)+nx+1:(i+1)*(nx+nu));

        mem2.gx(:,i+1) = gp(1:nx,i+1);
        mem2.gu(:,i+1) = gp(nx+1:end,i+1);
        
        mem2.Cgx(:,i*nx+1:(i+1)*nx)=Ccx(:,i*(nx+nu)+1:i*(nx+nu)+nx);
        mem2.Cgu(:,i*nu+1:(i+1)*nu)=Ccx(:,i*(nx+nu)+nx+1:(i+1)*(nx+nu));
        
        mem2.lb_dx(i*nbx+1:(i+1)*nbx,1) = mem.lb_dx(i*Npc*nbx+1:i*Npc*nbx+nbx);
        mem2.ub_dx(i*nbx+1:(i+1)*nbx,1) = mem.ub_dx(i*Npc*nbx+1:i*Npc*nbx+nbx);
    end
    mem2.Q(:,N2*nx+1:(N2+1)*nx) = mem.Q(:,N*nx+1:(N+1)*nx);
    mem2.gx(:,N2+1) = mem.gx(:,N+1);    
end


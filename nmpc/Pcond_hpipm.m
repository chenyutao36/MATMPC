function [mem2] = Pcond_hpipm(mem, settings, mem2, settings2)
    
    [Qp,Sp,Rp,qp,rp,Ap,Bp,bp,ubxp,lbxp,ubup,lbup,Cp,Dp,CNp,ugp,lgp,ds0p]=Partial_Condensing(mem,settings);

    mem2.gx=qp;
    mem2.gu=rp;
    mem2.A=Ap;
    mem2.B=Bp;
    mem2.a=bp;
    mem2.lb_dx=lbxp;
    mem2.ub_dx=ubxp;
    mem2.lb_du=lbup;
    mem2.ub_du=ubup;    
    mem2.Cgx=Cp;
    mem2.Cgu=Dp;
    mem2.CgN=CNp;
    mem2.lc=lgp;
    mem2.uc=ugp;
    mem2.ds0=ds0p;
                                
    for i=0:settings2.N-1
        mem2.Q(:,i*settings2.nx+1:(i+1)*settings2.nx) = Qp(:,i*settings2.nx+1:(i+1)*settings2.nx);
        
        mem2.R(:,i*settings2.nu+1:(i+1)*settings2.nu) = Rp(:,i*settings2.nu+1:(i+1)*settings2.nu);

    end
    mem2.Q(:,settings2.N*settings2.nx+1:(settings2.N+1)*settings2.nx) = Qp(:,settings2.N*settings2.nx+1:(settings2.N+1)*settings2.nx);
        
    for i=0:settings2.N-1
        mem2.S(:,i*settings2.nu+1:(i+1)*settings2.nu)=(Sp(:,i*settings2.nx+1:(i+1)*settings2.nx))';
    end
end


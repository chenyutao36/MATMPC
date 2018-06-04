function [mem2] = Pcond(mem, settings, mem2, settings2)

    nu = settings2.nu;
    
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
%         Qtmp = Qp(:,i*settings2.nx+1:(i+1)*settings2.nx);
%         Qtmp2 = Qtmp + Qtmp';
% %         Qtmp2(1:settings2.nx+1:end) = diag(Qtmp);
%         mem2.Q(:,i*settings2.nx+1:(i+1)*settings2.nx) = Qtmp2;
        mem2.Q(:,i*settings2.nx+1:(i+1)*settings2.nx) = Qp(:,i*settings2.nx+1:(i+1)*settings2.nx);
        
%         Rtmp = Rp(:,i*settings2.nu+1:(i+1)*settings2.nu);
% %         Rtmp = rot90(Rp(:,i*settings2.nu+1:(i+1)*settings2.nu), 2);
%         Rtmp2 = Rtmp + Rtmp';
%         Rtmp2(1:settings2.nu+1:end) = diag(Rtmp);
%         mem2.R(:,i*settings2.nu+1:(i+1)*settings2.nu) = Rtmp2;
        mem2.R(:,i*settings2.nu+1:(i+1)*settings2.nu) = Rp(:,i*settings2.nu+1:(i+1)*settings2.nu);

%         mem2.gu(i*nu+1:(i+1)*nu) = flipud(rp(i*nu+1:(i+1)*nu));
%         mem2.B(:,i*nu+1:(i+1)*nu) = fliplr(Bp(:,i*nu+1:(i+1)*nu));
%         mem2.lb_du(i*nu+1:(i+1)*nu) = flipud(lbup(i*nu+1:(i+1)*nu));
%         mem2.ub_du(i*nu+1:(i+1)*nu) = flipud(ubup(i*nu+1:(i+1)*nu));
%         mem2.Cgu(:,i*nu+1:(i+1)*nu) = rot90(Dp(:,i*nu+1:(i+1)*nu),2);
    end
%     Qtmp = Qp(:,settings2.N*settings2.nx+1:(settings2.N+1)*settings2.nx);
%     Qtmp2 = Qtmp + Qtmp';
%     Qtmp2(1:settings2.nx+1:end) = diag(Qtmp);
%     mem2.Q(:,settings2.N*settings2.nx+1:(settings2.N+1)*settings2.nx) = Qtmp2;
    mem2.Q(:,settings2.N*settings2.nx+1:(settings2.N+1)*settings2.nx) = Qp(:,settings2.N*settings2.nx+1:(settings2.N+1)*settings2.nx);
        
    for i=0:settings2.N-1
        mem2.S(:,i*settings2.nu+1:(i+1)*settings2.nu)=(Sp(:,i*settings2.nx+1:(i+1)*settings2.nx))';
%         mem2.S(:,i*settings2.nu+1:(i+1)*settings2.nu)=flipud((Sp(:,i*settings2.nx+1:(i+1)*settings2.nx))');
    end
end


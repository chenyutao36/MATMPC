function [mem2] = Pcond(mem, settings, mem2, settings2)
        
    [Hp,gp,Ccx,Ap,Bp,ap,lxc,uxc,Ccg,lgc,ugc]=partial_condensing_default(mem, settings);
    
    sparse2full(mem2, mem, settings2, settings, Hp, gp, Ccx, Ccg, lxc, uxc, lgc, ugc);
    
    mem2.A=Ap;
    mem2.B=Bp;
    mem2.a=ap;
    mem2.lb_du=mem.lb_du;
    mem2.ub_du=mem.ub_du;    
    mem2.ds0=mem.ds0;
    mem2.CgN=mem.CgN;
    mem2.Cx=mem.Cx;

end


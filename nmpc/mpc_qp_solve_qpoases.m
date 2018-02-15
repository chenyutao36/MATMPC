function [cpt_qp, mem] = mpc_qp_solve_qpoases(sizes,mem)

    nu=sizes.nu;
    N=sizes.N;   
        
    H=mem.Hc;
    g=mem.gc;
    C=mem.Cc;
    lu=mem.lb_du;
    uu=mem.ub_du;
    lc=mem.lcc;
    uc=mem.ucc;
    
    options = qpOASES_options('MPC');
%     options = qpOASES_options('default');
    if mem.warm_start==0
        
               
        [QP,solution,fval,exitflag,iterations,multiplier,auxOutput] = qpOASES_sequence('i',H,g,C,lu,uu,lc,uc, options); 
        mem.warm_start=QP;
    else
        QP=mem.warm_start;
        if mem.hot_start==0
            [solution,fval,exitflag,iterations,multiplier,auxOutput] = qpOASES_sequence('m',QP,H,g,C,lu,uu,lc,uc, options);
        end
        if mem.hot_start==1
           [solution,fval,exitflag,iterations,multiplier,auxOutput] = qpOASES_sequence('h',QP,g,lu,uu,lc,uc, options);
        end
    end
    mem.mu_u_new  = - multiplier(1:N*nu);
    mu_vec   = - multiplier(N*nu+1:end);
    cpt_qp   = auxOutput.cpuTime*1e3;
            
    Recover(mem, sizes, solution, mu_vec);
end

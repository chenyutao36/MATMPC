function [solution,mu_vec,cpt_qp,mem] = mpc_qp_solve_dense(H,g,C,lu,uu,lc,uc,sizes,mem)

    nu=sizes.nu;
    N=sizes.N;   
    options = qpOASES_options('MPC');
    if mem.qpoases.warm_start==0
        [QP,solution,fval,exitflag,iterations,multiplier,auxOutput] = qpOASES_sequence('i',H,g,C,lu,uu,lc,uc, options); 
        mem.qpoases.warm_start=QP;
    else
        QP=mem.qpoases.warm_start;
        if mem.qpoases.hot_start==0
            [solution,fval,exitflag,iterations,multiplier,auxOutput] = qpOASES_sequence('m',QP,H,g,C,lu,uu,lc,uc, options);
        end
        if mem.qpoases.hot_start==1
           [solution,fval,exitflag,iterations,multiplier,auxOutput] = qpOASES_sequence('h',QP,g,lu,uu,lc,uc, options);
        end
    end           
    mu_vec   = - multiplier(N*nu+1:end);
    cpt_qp   = auxOutput.cpuTime*1e3;
            

end

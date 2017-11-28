function [solution,mu_vec,cpt_qp,mem] = mpc_qp_solve_dense(H,g,dB,B,sizes, mem)

    nu=sizes.nu;   
    N=sizes.N;   

    
    if mem.qpoases.warm_start==0
        [QP,solution,fval,exitflag,iterations,multiplier,auxOutput] = qpOASES_sequence('i',H,g,dB,[],[],[],-B); 
        mem.qpoases.warm_start=QP;
    else
        QP=mem.qpoases.warm_start;
        if mem.hot_start==0
             [solution,fval,exitflag,iterations,multiplier,auxOutput] = qpOASES_sequence('m',QP,H,g,dB,[],[],[],-B);
        end
        if mem.hot_start==1
             [solution,fval,exitflag,iterations,multiplier,auxOutput] = qpOASES_sequence('h',QP,g,[],[],[],-B);
        end
    end           
    mu_vec   = - multiplier(N*nu+1:end);
    cpt_qp   = auxOutput.cpuTime*1e3;
            

end

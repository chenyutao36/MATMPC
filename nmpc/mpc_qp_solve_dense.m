function [solution,mu_vec,cpt_qp,input] = mpc_qp_solve_dense(H,g,dB,B,input,sizes)

    nu=sizes.nu;   
    N=sizes.N;   

    
    if input.opt.condensing_matrix.QP==0 
        [QP,solution,fval,exitflag,iterations,multiplier,auxOutput] = qpOASES_sequence('i',H,g,dB,[],[],[],-B); 
        input.opt.condensing_matrix.QP=QP;
    else
        QP=input.opt.condensing_matrix.QP;
        if strcmp(input.opt.hotstart,'no')
             [solution,fval,exitflag,iterations,multiplier,auxOutput] = qpOASES_sequence('m',QP,H,g,dB,[],[],[],-B);
        end
        if strcmp(input.opt.hotstart,'yes')
             [solution,fval,exitflag,iterations,multiplier,auxOutput] = qpOASES_sequence('h',QP,g,[],[],[],-B);
        end
    end           
    mu_vec   = - multiplier(N*nu+1:end);
    cpt_qp   = auxOutput.cpuTime*1e3;
            

end

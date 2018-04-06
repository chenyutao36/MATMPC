function [output, mem] = mpc_nmpcsolver(input, settings, mem, opt)

    tic;

    i=0;
    KKT = 1e8;
    
    CPT.SHOOT=0;
    CPT.COND=0;
    CPT.QP=0;
   
    while(i < mem.sqp_maxit  &&  KKT > mem.kkt_lim ) % RTI or multiple call
        
        %% ----------- QP Preparation
       
        tshoot = tic;
        qp_generation(input, settings, mem);
        tSHOOT = toc(tshoot)*1e3; 
        
        switch opt.condensing
            case 'default_full'
                tcond=tic;
                Condensing(mem, settings);
                tCOND=toc(tcond)*1e3;
            case 'hpipm_full'
                tcond=tic;
                condensing_hpipm(mem, settings);
                tCOND=toc(tcond)*1e3;
            case 'no'
                tCOND = 0;
        end
                                            
        %% ----------  Solving QP
        switch opt.qpsolver
            case 'qpoases'              
                [tQP,mem] = mpc_qp_solve_qpoases(settings,mem, opt);
            case 'qore'
                [tQP,mem] = mpc_qp_solve_qore(settings,mem, opt);
            case 'quadprog'
                [tQP,mem] = mpc_qp_solve_quadprog(settings,mem);
            case 'hpipm_sparse'               
                tqp=tic;
                hpipm_sparse(mem,settings);
                tQP = toc(tqp)*1e3;
            case 'hpipm_pcond'               
                tqp=tic;
                hpipm_pcond(mem,settings);
                tQP = toc(tqp)*1e3;
        end
        

        %% ---------- Line search

        Line_search(mem, input, settings);
                
        %% ---------- KKT calculation 
        
        [eq_res, ineq_res, KKT] = solution_info(input, settings, mem);
        
        %% ---------- Multiple call management and convergence check
                        
        CPT.SHOOT=CPT.SHOOT+tSHOOT;
        CPT.COND=CPT.COND+tCOND;
        CPT.QP=CPT.QP+tQP;
        
        i=i+1;
        
    end

    output.info.cpuTime=toc*1e3;   % Total CPU time for the current sampling instant
    
    output.x=input.x;
    output.u=input.u;   
    output.lambda=input.lambda;
    output.mu=input.mu;
    output.muN=input.muN;
    output.mu_u=input.mu_u;

    output.info.iteration_num=i;    
    output.info.kktValue=KKT;
    output.info.eq_res=eq_res;
    output.info.ineq_res=ineq_res;
    output.info.shootTime=CPT.SHOOT;
    output.info.condTime=CPT.COND;
    output.info.qpTime=CPT.QP;

end


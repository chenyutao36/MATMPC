function [output, mem] = mpc_nmpcsolver(input, settings, mem)

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
                      
        tcond=tic;
        Condensing(mem, settings);        
        tCOND=toc(tcond)*1e3;
                
%         %% ----------  Solving QP
        [tQP,mem] = mpc_qp_solve_qpoases(settings,mem);
        
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

    output.info.iteration_num=i;    
    output.info.kktValue=KKT;
    output.info.eq_res=eq_res;
    output.info.ineq_res=ineq_res;
    output.info.shootTime=CPT.SHOOT;
    output.info.condTime=CPT.COND;
    output.info.qpTime=CPT.QP;

end


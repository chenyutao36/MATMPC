function [input,mem,cpt,StopCrit,OBJ] = mpc_nmpcsolver_simulink(input, settings, mem, opt)

    mem.sqp_it=0;
    mem.alpha =1;
    mem.obj=0;
    StopCrit = 2*mem.kkt_lim;
    
    tic;
    while(mem.sqp_it < mem.sqp_maxit  &&  StopCrit > mem.kkt_lim && mem.alpha>1E-8 )
               
        if strcmp(opt.qpsolver, 'qpoases_mb')
            qp_generation_mb(input, settings, mem);
        else
            qp_generation(input, settings, mem);
        end

        switch opt.condensing
             case 'default_full'  
                 if ~strcmp(opt.qpsolver, 'qpoases_mb')
                    Condensing(mem, settings);
                 else
                    Condensing_mb(mem, settings);
                 end
             case 'blasfeo_full'               
                 Condensing_Blasfeo(mem, settings);               
             case 'no'

         end

        switch opt.qpsolver
            case 'qpoases'              
                [~,mem] = mpc_qp_solve_qpoases(settings,mem);
            case 'qpoases_mb'              
                [~,mem] = mpc_qp_solve_qpoases_mb(settings,mem, opt);
            case 'hpipm_sparse'               
                hpipm_sparse(mem,settings);
            case 'hpipm_pcond'                             
                hpipm_pcond(mem,settings);

        end

        Line_search(mem, input, settings);

        [eq_res, ineq_res, KKT, OBJ] = solution_info(input, settings, mem);

        StopCrit = max([eq_res, ineq_res, KKT]);

        mem.sqp_it=mem.sqp_it+1;
    
    end
    cpt=toc*1e3;
    
end


function [input,mem,cpt] = mpc_nmpcsolver_simulink(input, settings, mem, opt)

    tic;
               
    qp_generation(input, settings, mem);

    switch opt.condensing
         case 'default_full'           
             Condensing(mem, settings);
         case 'blasfeo_full'               
             Condensing_Blasfeo(mem, settings);               
         case 'no'
             
     end
           
    switch opt.qpsolver
        case 'qpoases'              
            [~,mem] = mpc_qp_solve_qpoases(settings,mem);
        case 'hpipm_sparse'               
            hpipm_sparse(mem,settings);
        case 'hpipm_pcond'                             
            hpipm_pcond(mem,settings);
                
    end
    
    Line_search(mem, input, settings);
    
%     [eq_res, ineq_res, KKT] = solution_info(input, settings, mem);

    cpt=toc*1e3;
    
end


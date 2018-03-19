function [input,mem,cpt] = mpc_nmpcsolver_simulink(input, settings, mem)

    tic;
               
    qp_generation(input, settings, mem);

    Condensing(mem, settings);
           
    [~,mem] = mpc_qp_solve_qpoases(settings,mem);
    
    Line_search(mem, input, settings);

    cpt=toc*1e3;
    
end


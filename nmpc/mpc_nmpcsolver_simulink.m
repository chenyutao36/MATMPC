function [input,mem,cpt,StopCrit,OBJ,lineSearchAccepted] = mpc_nmpcsolver_simulink(input, settings, mem, opt)

    legacyInput = ~isfield(input, 'z');
    if legacyInput
        if settings.nz ~= 0
            error('MATMPC:MissingAlgebraicState', 'The Simulink controller requires input.z for a DAE model.');
        end
        % The R2017b pendulum chart omits the empty algebraic state field.
        input.z = zeros(0, settings.N);
    end

    mem.sqp_it=0;
    mem.alpha =1;
    mem.obj=0;
    mem.nonuniform_grid = isfield(opt, 'nonuniform_grid') && logical(opt.nonuniform_grid);
    if mem.nonuniform_grid && strcmp(opt.qpsolver, 'qpoases_mb')
        error('MATMPC:IncompatibleGrid', 'Move blocking requires a uniform shooting grid.');
    end
    StopCrit = 2*mem.kkt_lim;
    
    tic;
    while(mem.sqp_it < mem.sqp_maxit  &&  StopCrit > mem.kkt_lim && mem.alpha>1E-8 )
               
        if strcmp(opt.qpsolver, 'qpoases_mb')
            qp_generation_mb(input, settings, mem);
        else
            mem = prepare_shooting_qp(input, settings, mem);
        end

        switch opt.condensing
             case 'default_full'  
                 if ~strcmp(opt.qpsolver, 'qpoases_mb')
                    mem = condense_default_qp(mem, settings);
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
            otherwise
                error('MATMPC:UnsupportedQPSolver', 'Unsupported Simulink QP solver "%s".', opt.qpsolver);
        end

        [input, mem, eq_res, ineq_res, KKT, OBJ] = mpc_sqp_step(input, settings, mem);

        StopCrit = max([eq_res, ineq_res, KKT]);

        mem.sqp_it=mem.sqp_it+1;
    
    end
    cpt=toc*1e3;
    lineSearchAccepted = mem.alpha > 0;

    if legacyInput
        input = rmfield(input, 'z');
    end
end


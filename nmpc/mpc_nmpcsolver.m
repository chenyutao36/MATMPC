function [output, mem] = mpc_nmpcsolver(input, settings, mem, opt)


    mem.sqp_it=0;
    mem.alpha =1;
    StopCrit = 2*mem.kkt_lim;
    
    CPT.SHOOT=0;
    CPT.COND=0;
    CPT.QP=0;

    tic;
  
    while(mem.sqp_it < mem.sqp_maxit  &&  StopCrit > mem.kkt_lim && mem.alpha>1E-4 ) % RTI or multiple call
        
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
            case 'blasfeo_full'
                tcond=tic;
                Condensing_Blasfeo(mem, settings);
                tCOND=toc(tcond)*1e3;
            case 'partial_condensing'
                tcond=tic;
                mem.mem2 = Pcond(mem, settings, mem.mem2, mem.settings2);
                tCOND=toc(tcond)*1e3;                
            case 'no'
                tCOND = 0;
        end
                                            
        %% ----------  Solving QP
        switch opt.qpsolver
            case 'qpoases'              
                [tQP,mem] = mpc_qp_solve_qpoases(settings,mem);
            case 'qore'
                [tQP,mem] = mpc_qp_solve_qore(settings,mem);
            case 'quadprog_dense'
                [tQP,mem] = mpc_qp_solve_quadprog(settings,mem);
            case 'hpipm_sparse'               
                tqp=tic;
                hpipm_sparse(mem,settings);
                tQP = toc(tqp)*1e3;
            case 'hpipm_pcond'               
                tqp=tic;
                hpipm_pcond(mem,settings);
                tQP = toc(tqp)*1e3;
            case 'ipopt_dense'               
                [tQP,mem] = mpc_qp_solve_ipopt_dense(settings,mem);
            case 'ipopt_sparse'               
                [tQP,mem] = mpc_qp_solve_ipopt_sparse(settings,mem);
                
            case 'ipopt_partial_sparse'
                [tQP, mem] = mpc_qp_solve_ipopt_partial_sparse(settings,mem.settings2,mem, mem.mem2);               
                
            case 'qpdunes'
                [tQP, mem] = mpc_qp_solve_qpdunes(settings,mem);
                
%                 mem2.iter=mem.iter;
%                 mem2.dunes.qpOptions=mem.dunes.qpOptions;
%                 nw = (settings2.N+1)*settings2.nx+settings2.N*settings2.nu;
%                 mem2.dunes.H=zeros(settings2.nx+settings2.nu,(settings2.nx+settings2.nu)*settings2.N);
%                 mem2.dunes.P=zeros(settings2.nx,settings2.nx);
%                 mem2.dunes.C=zeros(settings2.nx,(settings2.nx+settings2.nu)*settings2.N);
%                 mem2.dunes.c=zeros(settings2.nx, settings2.N);
%                 mem2.dunes.g=zeros(nw,1);
%                 mem2.dunes.zLow = -inf(nw,1);
%                 mem2.dunes.zUpp = inf(nw,1);
%                 mem2.dunes.D = zeros(settings2.nc,(settings2.nx+settings2.nu)*settings2.N+settings2.nx);
%                 mem2.dunes.dLow = zeros(settings2.nc*(settings2.N+1));
%                 mem2.dunes.dUpp = zeros(settings2.nc*(settings2.N+1));
%                 [tQP, mem2] = mpc_qp_solve_qpdunes(settings2,mem2);
        end
        

        %% ---------- Line search

        Line_search(mem, input, settings);
                
        %% ---------- KKT calculation 
        
        [eq_res, ineq_res, KKT] = solution_info(input, settings, mem);
        
        StopCrit = max([eq_res, ineq_res, KKT]);
        
        %% ---------- Multiple call management and convergence check
                        
        CPT.SHOOT=CPT.SHOOT+tSHOOT;
        CPT.COND=CPT.COND+tCOND;
        CPT.QP=CPT.QP+tQP;
        
        mem.sqp_it=mem.sqp_it+1;
              
    end

    output.info.cpuTime=toc*1e3;   % Total CPU time for the current sampling instant
    
    output.x=input.x;
    output.u=input.u;   
    output.lambda=input.lambda;
    output.mu=input.mu;
    output.mu_x=input.mu_x;
    output.mu_u=input.mu_u;

    output.info.iteration_num=mem.sqp_it;      
    output.info.kktValue=KKT;
    output.info.OptCrit = StopCrit;
    output.info.eq_res=eq_res;
    output.info.ineq_res=ineq_res;
    output.info.shootTime=CPT.SHOOT;
    output.info.condTime=CPT.COND;
    output.info.qpTime=CPT.QP;

end


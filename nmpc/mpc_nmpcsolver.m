function [output, mem] = mpc_nmpcsolver(input,settings, mem)
    tic;

%     i=0;
%     KKT = 1e8;
%     
%     CPT.SHOOT=0;
%     CPT.COND=0;
%     CPT.QP=0;
   
%     while(i < mem.sqp_maxit  &&  KKT > mem.kkt_lim ) % RTI or multiple call
        
        %% ----------- QP Preparation

        tshoot = tic;
        qp_generation(input, settings, mem);
%         qp_generation_cmon(input, settings, mem);
        tSHOOT = toc(tshoot)*1e3; 
        
               
        tcond=tic;
        Condensing(mem, settings);
        tCOND=toc(tcond)*1e3;
        
%         %% ----------  Solving QP
        [du, mu_vec,tQP,mem] = mpc_qp_solve_dense(settings,mem);
% 
        Recover(mem, settings, du, mu_vec);

        %% hpipm test
        
%         [dz_hpipm, dxN_hpipm, lambda_hpipm, mu_hpipm, muN_hpipm, tQP_hpipm] = mpc_qp_solve_sparse(settings,mem);
%         tCOND = 0;        
%         err = [norm(mem.dz-dz_hpipm), norm(mem.dxN-dxN_hpipm), norm(mem.lambda_new(:,2:end)-lambda_hpipm), norm(mem.mu_new-mu_hpipm) norm(mem.muN_new-muN_hpipm)]
        
%         dz = dz_hpipm;
%         dxN = dxN_hpipm;
%         lambda = lambda_hpipm;
%         mu = mu_hpipm;
%         muN = muN_hpipm;
%         tQP = tQP_hpipm;

        %% ---------- Line search

        Line_search(mem, input, settings);

        %% ---------- KKT calculation 
        
        [eq_res, ineq_res, KKT] = solution_info(input, settings, mem);
%         eq_res = 0; ineq_res = 0; KKT = 0;
        
        %% ---------- Multiple call management and convergence check
                        
%         CPT.SHOOT=CPT.SHOOT+tSHOOT;
%         CPT.COND=CPT.COND+tCOND;
%         CPT.QP=CPT.QP+tQP;
    
        CPT.SHOOT=tSHOOT;
        CPT.COND=tCOND;
        CPT.QP=tQP;
               
%         i=i+1;
        
%         if norm(mem.dz,1)+norm(mem.dxN,1)<1e-4
%             break;
%         end
%     end

    output.info.cpuTime=toc*1e3;   % Total CPU time for the current sampling instant
    
    output.z=input.z;
    output.xN=input.xN;   
    output.lambda=input.lambda;
    output.mu=input.mu;
    output.muN=input.muN;

%     output.info.iteration_num=i;    
    output.info.kktValue=KKT;
    output.info.eq_res=eq_res;
    output.info.ineq_res=ineq_res;
    output.info.shootTime=CPT.SHOOT;
    output.info.condTime=CPT.COND;
    output.info.qpTime=CPT.QP;

end


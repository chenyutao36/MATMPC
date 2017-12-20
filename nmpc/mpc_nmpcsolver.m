function [output, mem] = mpc_nmpcsolver(input,settings, mem)
    tic;

%     i=0;
%     KKT = 1e8;
    
%     CPT.INT=0;
%     CPT.SENS=0;
%     CPT.COND=0;
%     CPT.QP=0;
   
%     while(i < mpc_callnum  &&  KKT > kkt_lim) % RTI or multiple call
        
        %% ----------- QP Preparation

        tshoot = tic;
        qp_generation(input, settings, mem);
        tSHOOT = toc(tshoot)*1e3; 
        
               
        tcond=tic;
        Condensing(mem, settings);
        tCOND=toc(tcond)*1e3;
        
%         %% ----------  Solving QP
        [du,mu_vec,tQP,mem] = mpc_qp_solve_dense(settings,mem);
% 
        Recover(du,mu_vec,mem,settings);

        %% hpipm test
        
%         [dz_hpipm, dxN_hpipm, lambda_hpipm, mu_hpipm, muN_hpipm, tQP_hpipm] = mpc_qp_solve_sparse(Q_h,S,R,A,B,Cx,Cu,gx,gu,a,ds0,lc,uc,lb_du,ub_du,CxN,settings,mem);
%         tCOND = 0;        
%         err = [norm(dz-dz_hpipm), norm(lambda(:,2:end)-lambda_hpipm), norm(mu-mu_hpipm)]
        
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
        
        
        %% ---------- Multiple call management and convergence check
                        
%         CPT.INT=CPT.INT+tINT;
%         CPT.SENS=CPT.SENS+tSENS;
%         CPT.COND=CPT.COND+tcond;
%         CPT.QP=CPT.QP+info.cpt_qp;
    
        CPT.SHOOT=tSHOOT;
        CPT.COND=tCOND;
        CPT.QP=tQP;
               
%         i=i+1;
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


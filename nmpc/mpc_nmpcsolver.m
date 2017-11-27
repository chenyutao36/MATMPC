function [output] = mpc_nmpcsolver(input,settings, sim_erk_mem, sim_irk_mem)
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
%         [Q_h,S,R,A,B,Cx,Cu,gx,gu,c,a,ds0] = qp_generation_mex(input, settings);
        [Q_h,S,R,A,B,Cx,Cu,gx,gu,c,a,ds0] = qp_generation_mex_erk(input, settings, sim_erk_mem);
%         [Q_h,S,R,A,B,Cx,Cu,gx,gu,c,a,ds0] = qp_generation_mex_irk(input, settings, sim_irk_mem);
        tSHOOT = toc(tshoot)*1000; 
        
        tcond=tic;
        [Hc,gc, Cc,cc]=Condensing_mex_new(A,B,Q_h,S,R,Cx,Cu,ds0,a,c,gx,gu,settings);
        tCOND=toc(tcond)*1e3;
        
        %% ----------  Solving QP
        [du,mu_vec,tQP,input] = mpc_qp_solve_dense(Hc,gc,Cc,cc,input,settings);

        [dz, dxN, lambda, mu, muN] = Recover_mex(Q_h,S,A,B,Cx,a,gx,du,ds0,mu_vec,settings);

        %% ---------- Line search

        alpha = 1;
        [z,xN,lambda,mu,muN] = Line_search_mex(dz, dxN, lambda, mu, muN, alpha, input, settings);

        %% ---------- KKT calculation 
        
        [eq_res, ineq_res, KKT] = solution_info(lambda, mu, muN, ds0, input, settings);
        %% ---------- Multiple call management and convergence check
                        
%         CPT.INT=CPT.INT+tINT;
%         CPT.SENS=CPT.SENS+tSENS;
%         CPT.COND=CPT.COND+tcond;
%         CPT.QP=CPT.QP+info.cpt_qp;
        
        CPT.SHOOT=tSHOOT;
        CPT.COND=tCOND;
        CPT.QP=tQP;
               
%         input.opt.meritfun.mu_merit = mu_merit;
        
%         i=i+1;
%     end

    output.info.cpuTime=toc*1e3;   % Total CPU time for the current sampling instant
    
    output.z=z;
    output.xN=xN;   
    output.lambda=lambda;
    output.mu=mu;
    output.muN=muN;

    % output.info.iteration_num=i;    
    output.info.kktValue=KKT;
    output.info.eq_res=eq_res;
    output.info.ineq_res=ineq_res;
    output.info.shootTime=CPT.SHOOT;
    output.info.condTime=CPT.COND;
    output.info.qpTime=CPT.QP;
    
    if strcmp(input.opt.qpsolver,'qpoases')
        output.opt.condensing_matrix.QP=input.opt.condensing_matrix.QP;
    end

%     output.meritfun=input.opt.meritfun;
end


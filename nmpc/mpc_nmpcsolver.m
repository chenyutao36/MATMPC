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
        [Q_h,S,R,A,B,Cx,Cu,gx,gu,a,ds0,lc,uc,lb_du,ub_du,CxN] = qp_generation(input, settings, mem);
        tSHOOT = toc(tshoot)*1e3; 
        
        
        
        tcond=tic;
        [Hc,gc, Cc, lcc, ucc] = Condensing(A,B,Q_h,S,R,Cx,Cu,CxN,ds0,a,gx,gu,lc,uc,settings);
        tCOND=toc(tcond)*1e3;
        
        %% ----------  Solving QP
        [du,mu_vec,tQP,mem] = mpc_qp_solve_dense(Hc,gc,Cc,lb_du,ub_du,lcc,ucc,settings,mem);

        [dz, dxN, lambda, mu, muN] = Recover(Q_h,S,A,B,Cx,CxN,a,gx,du,ds0,mu_vec,settings);

%         dz = dz_hpipm;
%         dxN = dxN_hpipm;
%         lambda = lambda_hpipm;
%         mu = mu_hpipm;
%         muN = muN_hpipm;

        %% hpipm test
        
        [dz_hpipm, dxN_hpipm, lambda_hpipm, mu_hpipm, muN_hpipm, tQP_hpipm] = mpc_qp_solve_sparse(Q_h,S,R,A,B,Cx,Cu,gx,gu,a,ds0,lc,uc,lb_du,ub_du,CxN,settings,mem);
        tCOND = 0;        
        err = [norm(dz-dz_hpipm), norm(lambda(:,2:end)-lambda_hpipm), norm(mu-mu_hpipm)]

        %% ---------- Line search

        [z,xN,lambda,mu,muN] = Line_search(dz, dxN, lambda, mu, muN, input, settings);

        %% ---------- KKT calculation 
        
        [eq_res, ineq_res, KKT] = solution_info(lambda, mu, muN, ds0, input, settings, mem);
        
        
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

%     output.info.iteration_num=i;    
    output.info.kktValue=KKT;
    output.info.eq_res=eq_res;
    output.info.ineq_res=ineq_res;
    output.info.shootTime=CPT.SHOOT;
    output.info.condTime=CPT.COND;
    output.info.qpTime=CPT.QP;

end


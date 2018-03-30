function [cpt_qp, mem] = mpc_qp_solve_qore(sizes,mem)

    nu=sizes.nu;
    nc=sizes.nc;
    ncN=sizes.ncN;
    N=sizes.N;   
           
    if ~isfield(mem,'qore_id')              
        [err, mem.qore_id] = QPDenseNew(N*nu, N*nc+ncN);
        QPDenseSetInt(mem.qore_id, 'prtfreq', -1);
        QPDenseSetData(mem.qore_id, (mem.Cc)', mem.Hc);
        
        t1 = tic;
        QPDenseOptimize(mem.qore_id, [mem.lb_du;mem.lcc], [mem.ub_du;mem.ucc], mem.gc);
        cpt_qp = toc(t1)*1e3;
    else
        QPDenseUpdateMatrices(mem.qore_id, (mem.Cc)', mem.Hc);
        
        t1 = tic;
        QPDenseOptimize(mem.qore_id, [mem.lb_du;mem.lcc], [mem.ub_du;mem.ucc], mem.gc);
        cpt_qp = toc(t1)*1e3;
    end
        
    pri_sol = QPDenseGetDblVector(mem.qore_id, 'primalsol');
    dual_sol = QPDenseGetDblVector(mem.qore_id, 'dualsol');
    
    mem.du = pri_sol(1:N*nu);
    mem.mu_u_new = dual_sol(1:N*nu);
    mu_vec = dual_sol(N*nu+1:end);
            
    Recover(mem, sizes, mu_vec);
end

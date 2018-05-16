function [cpt_qp, mem] = mpc_qp_solve_qore(sizes,mem)

    nu=sizes.nu;
    nc=sizes.nc;
    ncN=sizes.ncN;
    nbx=sizes.nbx;
    N=sizes.N; 
               
    if mem.qore_id==-1            
        [err, mem.qore_id] = QPDenseNew(N*nu, N*nbx+N*nc+ncN);
        QPDenseSetInt(mem.qore_id, 'prtfreq', -1);
        QPDenseSetData(mem.qore_id, [mem.Ccx;mem.Ccg]', mem.Hc);
        
        t1 = tic;
        QPDenseOptimize(mem.qore_id, [mem.lb_du;mem.lxc;mem.lcc], [mem.ub_du;mem.uxc;mem.ucc], mem.gc);
        cpt_qp = toc(t1)*1e3;
    else
        QPDenseUpdateMatrices(mem.qore_id, [mem.Ccx;mem.Ccg]', mem.Hc);
        
        t1 = tic;
        QPDenseOptimize(mem.qore_id, [mem.lb_du;mem.lxc;mem.lcc], [mem.ub_du;mem.uxc;mem.ucc], mem.gc);
        cpt_qp = toc(t1)*1e3;
    end
        
    pri_sol = QPDenseGetDblVector(mem.qore_id, 'primalsol');
    dual_sol = QPDenseGetDblVector(mem.qore_id, 'dualsol');
    
    mem.du = reshape(pri_sol(1:N*nu),[nu N]);
    mem.mu_u_new = -dual_sol(1:N*nu);
    mem.mu_x_new = -dual_sol(N*nu+1:N*nu+N*nbx);
    mem.mu_new = -dual_sol(N*nu+N*nbx+1:end);
                
    Recover(mem, sizes);
end

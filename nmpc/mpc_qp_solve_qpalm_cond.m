function [cpt_qp, mem] = mpc_qp_solve_qpalm_cond(sizes,mem)
    
    nu=sizes.nu;
    N=sizes.N;  
    nbx=sizes.nbx;
                   
%     if mem.qpalm_settings.warm_start==0               
%         mem.qpalm_solver.setup(mem.Hc,mem.gc, [eye(N*nu);mem.Ccx;mem.Ccg], [mem.lb_du;mem.lxc;mem.lcc],[mem.ub_du;mem.uxc;mem.ucc], mem.qpalm_settings); 
%         sol = mem.qpalm_solver.solve();
%         mem.qpalm_settings.warm_start=1;   
%     else
%         mem.qpalm_solver.setup(mem.Hc,mem.gc, [eye(N*nu);mem.Ccx;mem.Ccg], [mem.lb_du;mem.lxc;mem.lcc],[mem.ub_du;mem.uxc;mem.ucc], ...
%             mem.du, [mem.mu_u_new;mem.mu_x_new;mem.mu_new], mem.qpalm_settings); 
%         sol = mem.qpalm_solver.solve();
%     end
    mem.qpalm_solver=qpalm;
    mem.qpalm_settings=mem.qpalm_solver.default_settings();

    mem.qpalm_solver.setup(mem.Hc,mem.gc, [eye(N*nu);mem.Ccx;mem.Ccg], [mem.lb_du;mem.lxc;mem.lcc],[mem.ub_du;mem.uxc;mem.ucc], ...
            mem.du, [mem.mu_u_new;mem.mu_x_new;mem.mu_new], mem.qpalm_settings); 
    S = 'mem.qpalm_solver.solve()';
    [T,sol] = evalc(S);
    
    assert(strcmp(sol.info.status,'solved'), ['QP Error: ' sol.info.status]);
                
    mem.du(:) = reshape(sol.x, [nu,N]);
    mem.mu_u_new(:)  = sol.y(1:N*nu); 
    mem.mu_x_new(:) = sol.y(N*nu+1:N*nu+N*nbx);
    mem.mu_new(:)   = sol.y(N*nu+N*nbx+1:end);
    cpt_qp   = sol.info.run_time*1e3;
                
    Recover(mem, sizes);
    
    mem.qpalm_solver.delete();
end

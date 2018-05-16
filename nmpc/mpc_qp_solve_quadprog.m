function [cpt_qp, mem] = mpc_qp_solve_quadprog(sizes,mem)

    nc=sizes.nc;
    ncN=sizes.ncN;
    N=sizes.N; 
    nbx=sizes.nbx;
    
    warning off;
    
    tqp = tic;
    [mem.du,fval,exitflag,output,lambda] = quadprog(mem.Hc,mem.gc,[mem.Ccx;mem.Ccg;-mem.Ccx;-mem.Ccg],...
                                                    [mem.uxc;mem.ucc;-mem.lxc;-mem.lcc],[],[],mem.lb_du, mem.ub_du,[],mem.quadprog_opt);
                                                
    mem.mu_u_new = lambda.upper - lambda.lower;
    mem.mu_x_new = lambda.ineqlin(1:N*nbx) - lambda.ineqlin(N*nbx+N*nc+ncN+1:N*nbx+N*nc+ncN+N*nbx);   
    mem.mu_new = lambda.ineqlin(N*nbx+1:N*nbx+N*nc+ncN) - lambda.ineqlin(N*nbx+N*nc+ncN+N*nbx+1:end);
    cpt_qp   = toc(tqp)*1e3;
            
    Recover(mem, sizes);
end

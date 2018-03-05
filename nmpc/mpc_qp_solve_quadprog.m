function [cpt_qp, mem] = mpc_qp_solve_quadprog(sizes,mem)

    nc=sizes.nc;
    ncN=sizes.ncN;
    N=sizes.N;   
    
    tqp = tic;
    [mem.du,fval,exitflag,output,lambda] = quadprog(mem.Hc,mem.gc,[mem.Cc;-mem.Cc],...
                                                    [mem.ucc;-mem.lcc],[],[],mem.lb_du, mem.ub_du,[],mem.quadprog_opt);
    mem.mu_u_new = lambda.upper - lambda.lower;       
    
    mu_vec   = lambda.ineqlin(1:N*nc+ncN) - lambda.ineqlin(N*nc+ncN+1:end);
    cpt_qp   = toc(tqp)*1e3;
            
    Recover(mem, sizes, mu_vec);
end

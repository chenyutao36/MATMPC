function [cpt_qp, mem] = mpc_qp_solve_hpipm_dense(settings,mem)


    tqp=tic;
    [mu_vec]=hpipm_dense(mem,settings);
    cpt_qp = toc(tqp)*1e3;

    Recover(mem, settings, mu_vec);
end




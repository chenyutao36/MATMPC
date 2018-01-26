%%
% dz and dxN are primal solutions, lambda_new and mu_new are dual solutions
% for equality and inequality constraints, obtained from qpoases-3.2.1,
% which is correct

% the mex interface is compiled using 
% mex GCC='/usr/bin/gcc-4.9' HPIPM_d_solve_ipm2_hard_ocp_qp.c /home/chen/Documents/Packages/hpipm/lib/libhpipm.a /home/chen/Documents/Packages/blasfeo/lib/libblasfeo.a -I/home/chen/Documents/Packages/hpipm/include -I/home/chen/Documents/Packages/blasfeo/include 
% the source file may slightly different from yours since I did some polish
% to it
%%

load myData;

[dz_hpipm, dxN_hpipm, lambda_hpipm, mu_hpipm, muN_hpipm, tQP_hpipm] = mpc_qp_solve_sparse(settings,mem);
      
err = [norm(mem.dz-dz_hpipm), norm(mem.dxN-dxN_hpipm), norm(mem.lambda_new(:,2:end)-lambda_hpipm), 
       norm(mem.mu_new-mu_hpipm) norm(mem.muN_new-muN_hpipm)]
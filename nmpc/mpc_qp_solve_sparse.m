function [dz, dxN, lambda, mu, muN, cpt_qp] = mpc_qp_solve_sparse(sizes,mem)

Q_h=mem.Q_h;
S=mem.S;
R=mem.R;
A=mem.A_sens;
B=mem.B_sens;
Cx=mem.Cx;
CxN=mem.CxN;
Cu=mem.Cu;
gx=mem.gx;
gu=mem.gu;
a=mem.a;
ds0=mem.ds0;
lc=mem.lc;
uc=mem.uc;
lb_du=mem.lb_du;
ub_du=mem.ub_du;


nx = sizes.nx;
nu = sizes.nu;
N = sizes.N;

nb = nu; % number of two-sided box constraints
ng  = sizes.nc;            % number of two-sided general constraints         
ngN = sizes.ncN; % number of two-sided general constraints on last stage

time_invariant = 0; % time_invariant (0) vs time_variant (1) problems

nbu = min(nb, nu);

AA = A;
BB = B;
bb = a;

% cost function

QQ = Q_h(:,1:N*nx);
Qf = Q_h(:,N*nx:end);
SS = S;
RR = R;
qq = gx(:,1:N);
qf = gx(:,N+1);
rr = gu;

% general constraints

CN = CxN;

CC = Cx;
DD = Cu;

% box constraints
llb = [reshape(lb_du,[nu,N]), lb_du(1:nu)];
uub = [reshape(ub_du,[nu,N]), ub_du(1:nu)];
% general constraints
llg = reshape(lc(1:N*ng),[ng,N]);
uug = reshape(uc(1:N*ng),[ng,N]);
% general constraints last stage
lgN = lc(N*ng+1:end,1);
ugN = uc(N*ng+1:end,1);;

% initial guess for states and inputs
x = zeros(nx, N+1); x(:,1) = ds0; % initial condition
u = zeros(nu, N);
mult_pi = zeros(nx,N);
mult_lam = zeros(2*(nb+ng)*N+2*(nb+ngN),1);

free_x0 = 0; % consider x0 as optimization variable
warm_start = 0; % read initial guess from x and u

mu0 = 2;        % max element in cost function as estimate of max multiplier
sb = -1;		% actual number of performed iterations
k_max = 20;		% maximim number of iterations
tol = 1e-8;		% tolerance in the duality measure
infos = zeros(5, k_max);
inf_norm_res = zeros(1, 4);

tic;
% [x,u,mult_pi,mult_lam,infos,inf_norm_res,kk] = HPIPM_d_solve_ipm2_hard_ocp_qp(sb, k_max, mu0, tol, N, nx, nu, nb, ng, ngN, time_invariant, free_x0, warm_start, AA, BB, bb, QQ, Qf, RR, SS, qq, qf, rr, llb, uub, CC, DD, llg, uug, CN, lgN, ugN, ds0);
kk = HPIPM_d_solve_ipm2_hard_ocp_qp(sb, k_max, mu0, tol, N, nx, nu, nb, ng, ngN, time_invariant, free_x0, warm_start, AA, BB, bb, QQ, Qf, RR, SS, qq, qf, rr, llb, uub, CC, DD, llg, uug, CN, lgN, ugN, x,u,mult_pi,mult_lam,infos,inf_norm_res);
cpt_qp = toc*1e3;

% fprintf('\n  alpha_aff     mu_aff       sigma        alpha        mu\n');
% infos(:,1:kk)'
% inf_norm_res

dz = [x(:,1:N);u];
dxN = x(:,N+1);

lambda = mult_pi;

tmp = 2*(N+1)*nu;
mu = reshape(mult_lam(tmp+1:tmp+N*ng),[ng, N]) - reshape(mult_lam(tmp+N*ng+ngN+1:tmp+N*ng+ngN+N*ng),[ng, N]);

muN = mult_lam(tmp+N*ng+1:tmp+N*ng+ngN) - mult_lam(tmp+N*ng+ngN+N*ng:tmp+N*ng+ngN+N*ng+ngN);
end
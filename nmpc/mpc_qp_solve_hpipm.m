function [cpt_qp, mem] = mpc_qp_solve_hpipm(settings,mem)

% mex HPIPM_d_solve_ipm2_hard_ocp_qp.c /home/chen/Documents/Packages/hpipm/lib/libhpipm.a /home/chen/Documents/Packages/blasfeo/lib/libblasfeo.a -I/home/chen/Documents/Packages/hpipm/include -I/home/chen/Documents/Packages/blasfeo/include  % linear algebra in BLASFEO

    nx = settings.nx;			% number of states
    nu = settings.nu;		    % number of inputs (controls)
    N = settings.N;				% horizon length

    nb  = settings.nbu;		% number of two-sided box constraints
    ng  = settings.nc;      % number of two-sided general constraints
    ngN = settings.ncN;     % number of two-sided general constraints on last stage

    time_invariant = 0; % time_invariant (0) vs time_variant (1) problems

    AA = mem.A;
    BB = mem.B;
    bb = mem.a;

    QQ = mem.Q(:,1:N*nx);
    Qf = mem.Q(:,N*nx+1:(N+1)*nx);
    SS = mem.S;
    RR = mem.R;
    qq = mem.gx(:,1:N);
    qf = mem.gx(:,N+1);
    rr = mem.gu;

    CC = mem.Cx;
    DD = mem.Cu;
    CN = mem.CN;

    llb = [reshape(mem.lb_du,nu,N),zeros(nu,1)];
    uub = [reshape(mem.ub_du,nu,N),zeros(nu,1)];

    llg =  reshape(mem.lc(1:N*ng),ng,N);
    uug =  reshape(mem.uc(1:N*ng),ng,N);

    % general constraints last stage
    lgN = mem.lc(N*ng+1:end);
    ugN = mem.uc(N*ng+1:end);

    % initial guess for states and inputs
    x = zeros(nx, N+1); x(:,1) = mem.ds0; % initial condition
    u = zeros(nu, N);
    mult_pi = zeros(nx,N+1);
    mult_lam = zeros(2*(nb+ng)*N+2*(nb+ngN),1);

    free_x0 = 0; % consider x0 as optimization variable
    warm_start = 0; % read initial guess from x and u

    mu0 = 2;        % max element in cost function as estimate of max multiplier
    kk = -1;		% actual number of performed iterations
    k_max = 20;		% maximim number of iterations
    tol = 1e-8;		% tolerance in the duality measure
    infos = zeros(5, k_max);
    inf_norm_res = zeros(1, 4);

    tqp=tic;
%     HPIPM_d_solve_ipm2_hard_ocp_qp(kk, k_max, mu0, tol, N, nx, nu, nb, ng, ngN, time_invariant, free_x0, warm_start, AA, BB, bb, QQ, Qf, RR, SS, qq, qf, rr, llb, uub, CC, DD, llg, uug, CN, lgN, ugN, x, u, mult_pi, mult_lam, inf_norm_res, infos);
    HPIPM_d_solve_ipm2_hard_ocp_qp(kk, k_max, mu0, tol, N, nx, nu, nb, ng, ngN, time_invariant, free_x0, warm_start, AA, BB, bb, QQ, Qf, RR, SS, qq, qf, rr, llb, uub, CC, DD, llg, uug, CN, lgN, ugN, x, u, mult_pi, mem.mu_u_new, mem.mu_new, mem.muN_new, inf_norm_res, infos, mult_lam);
    cpt_qp = toc(tqp)*1e3;

    mem.du = u;
    mem.dx = x;
    mem.lambda_new = mult_pi;
    
end




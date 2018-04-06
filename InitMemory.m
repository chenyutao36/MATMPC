function [mem] = InitMemory(settings, opt, input)

    Ts  = settings.Ts;       % Sampling time
    Ts_st = settings.Ts_st;  % Shooting interval
    s = settings.s;      % number of integration steps per interval
    nx = settings.nx;    % No. of states
    nu = settings.nu;    % No. of controls
    ny = settings.ny;    % No. of outputs (references)    
    nyN= settings.nyN;   % No. of outputs at terminal stage 
    np = settings.np;    % No. of parameters (on-line data)
    nc = settings.nc;    % No. of constraints
    ncN = settings.ncN;  % No. of constraints at terminal stage
    N     = settings.N;             % No. of shooting points

    %% memory
    mem = struct;
    mem.warm_start=0;
    mem.hot_start=0;
    if strcmp(opt.hotstart, 'yes')
        mem.hot_start=1;
    end
    
    switch opt.qpsolver
        case 'qpoases'   
            mem.qpsolver = 0;
            mem.qpoases_opt = qpOASES_options('MPC');
        case 'qore'
            addpath(genpath('/home/chen/Documents/Packages/QORE/Matlab_dense'));
            mem.qpsolver = 1;
        case 'quadprog'
            mem.quadprog_opt.Algorithm = 'interior-point-convex';
            mem.quadprog_opt.Display = 'off';
            mem.quadprog_opt.OptimalityTolerance = 1e-6;
            mem.quadprog_opt.ConstraintTolerance = 1e-6;
            mem.quadprog_opt.StepTolerance = 1e-6;
        case 'hpipm_sparse'
            mem.mu0=1e2;
            mem.max_qp_it = 100;
            mem.pred_corr = 1;
            mem.cond_pred_corr = 1;
        case 'hpipm_pcond'
            mem.mu0=1e2;
            mem.max_qp_it = 100;
            mem.pred_corr = 1;
            mem.cond_pred_corr = 1;
    end
          
    switch opt.integrator
        case 'ERK4-CASADI'
            mem.sim_method = 0;
        case 'ERK4'
            mem.sim_method = 1;
            mem.A_tab=[0, 0, 0, 0;
                       0.5, 0, 0, 0;
                       0, 0.5, 0, 0;
                       0, 0, 1, 0];
            mem.B_tab=[1/6, 1/3, 1/3, 1/6];
            mem.num_steps = s;
            mem.num_stages = 4;
            mem.h=Ts_st/mem.num_steps;
            mem.nx = nx;
            mem.nu = nu;
            mem.Sx = eye(nx);
            mem.Su = zeros(nx,nu);
        case 'IRK3'
            mem.sim_method = 2;
            mem.A_tab=[5/36,             2/9-sqrt(15)/15, 5/36-sqrt(15)/30;
                       5/36+sqrt(15)/24, 2/9            , 5/36-sqrt(15)/24;
                       5/36+sqrt(15)/30, 2/9+sqrt(15)/15, 5/36];
            mem.B_tab=[5/18;4/9;5/18];
            mem.num_steps = s;
            mem.num_stages = 3;
            mem.h= Ts_st/mem.num_steps;
            mem.nx = nx;
            mem.nu = nu;
            mem.Sx = eye(nx);
            mem.Su = zeros(nx,nu);
            mem.newton_iter = 8;
            mem.JFK = mem.h*[mem.B_tab(1)*eye(nx,nx), mem.B_tab(2)*eye(nx,nx), mem.B_tab(3)*eye(nx,nx)];
        otherwise 
            error('Please choose a correct integrator');       
    end
    
    mem.sqp_maxit = 1;           % maximum number of iterations for each sampling instant (for RTI, this is ONE)
    mem.kkt_lim = 1e-1;          % tolerance on optimality
    mem.mu_merit=0;              % initialize the parameter
    mem.eta=1e-4;                % merit function parameter
    mem.tau=0.8;                 % step length damping factor
    mem.mu_safty=1.1;            % constraint weight update factor (for merit function)
    mem.rho=0.5;                 % merit function parameter
    
    mem.A = zeros(nx,nx*N);
    mem.B = zeros(nx,nu*N);
    mem.Cx = zeros(nc,nx*N);
    mem.Cu = zeros(nc,nu*N);
    mem.CN = zeros(ncN,nx);
    mem.gx = zeros(nx,N+1);
    mem.gu = zeros(nu,N);
    mem.a = zeros(nx,N);
    mem.ds0 = zeros(nx,1);
    mem.lc = zeros(N*nc+ncN,1);
    mem.uc = zeros(N*nc+ncN,1);
    mem.lb_du = zeros(N*nu,1);
    mem.ub_du = zeros(N*nu,1);
    
    mem.Hc = zeros(N*nu,N*nu);
    mem.Cc = zeros(N*nc+ncN,N*nu);
    mem.gc = zeros(N*nu,1);
    mem.lcc = zeros(N*nc+ncN,1);
    mem.ucc = zeros(N*nc+ncN,1);
    
    mem.Cc_qore = zeros(N*nu,N*nc+ncN);

    mem.dx = zeros(nx,N+1);
    mem.du = zeros(nu,N);
    mem.lambda_new = zeros(nx,N+1);
    mem.mu_new = zeros(nc,N);
    mem.muN_new = zeros(ncN,1);
    mem.mu_u_new = zeros(N*nu,1);
    
    if strcmp(opt.lin_obj,'yes')
        mem.lin_obj = 1;
        
        [Jx, Ju] = Ji_fun('Ji_fun',zeros(nx,1),zeros(nu,1),zeros(np,1),zeros(ny,1), input.W);
        Qi = full(Jx'*Jx);
        for i=1:nx
            if Qi(i,i)<1e-8
                Qi(i,i)=1e-8;
            end
        end
        Si = full(Jx'*Ju);
        Ri = full(Ju'*Ju);
        mem.Q = repmat(Qi,1,N+1);
        mem.S = repmat(Si,1,N);
        mem.R = repmat(Ri,1,N);
        
        JN = JN_fun('JN_fun',zeros(nx,1),zeros(np,1),zeros(nyN,1), input.WN);
        Qf = full(JN'*JN);
        for i=1:nx
            if Qf(i,i)<1e-8
                Qf(i,i)=1e-8;
            end
        end
        mem.Q(:,N*nx+1:end) = Qf;
    else
        mem.lin_obj = 0;
        mem.Q = zeros(nx,nx*(N+1));
        mem.S = zeros(nx,nu*N);
        mem.R = zeros(nu,nu*N);
    end
    mem.reg = 1e-8;
              
end


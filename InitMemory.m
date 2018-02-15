function [input, mem] = InitMemory(settings, opt, input)

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
    if strcmp(opt.qpsolver,'qpoases')
        mem.qpoases.warm_start=0;
        mem.qpoases.hot_start=0;
        if strcmp(opt.hotstart, 'yes')
            mem.qpoases.hot_start=1;
        end
    end
      
    switch opt.integrator
        case 'ERK4-CASADI'
            mem.sim_method = 0;
        case 'ERK4'
            mem.sim_method = 1;
            mem.A=[0, 0, 0, 0;
                           0.5, 0, 0, 0;
                           0, 0.5, 0, 0;
                           0, 0, 1, 0];
            mem.B=[1/6, 1/3, 1/3, 1/6];
            mem.num_steps = s;
            mem.num_stages = 4;
            mem.h=Ts_st/mem.num_steps;
            mem.nx = nx;
            mem.nu = nu;
            mem.Sx = eye(nx);
            mem.Su = zeros(nx,nu);
        case 'IRK3'
            mem.sim_method = 2;
            mem.A=[5/36,             2/9-sqrt(15)/15, 5/36-sqrt(15)/30;
                   5/36+sqrt(15)/24, 2/9            , 5/36-sqrt(15)/24;
                   5/36+sqrt(15)/30, 2/9+sqrt(15)/15, 5/36];
            mem.B=[5/18;4/9;5/18];
            mem.num_steps = s;
            mem.num_stages = 3;
            mem.h= Ts_st/mem.num_steps;
            mem.nx = nx;
            mem.nu = nu;
            mem.Sx = eye(nx);
            mem.Su = zeros(nx,nu);
            mem.newton_iter = 8;
            mem.JFK = mem.h*[mem.B(1)*eye(nx,nx), mem.B(2)*eye(nx,nx), mem.B(3)*eye(nx,nx)];
        otherwise 
            error('Please choose a correct integrator');       
    end
    
    mem.sqp_maxit = 1;           % maximum number of iterations for each sampling instant (for RTI, this is ONE)
    mem.kkt_lim = 1;             % tolerance on optimality
    mem.mu_merit=0;              % initialize the parameter
    mem.eta=1e-4;                % merit function parameter
    mem.tau=0.8;                 % step length damping factor
    mem.mu_safty=1.1;            % constraint weight update factor (for merit function)
    mem.rho=0.5;                 % merit function parameter
    
    mem.A_sens = zeros(nx,nx*N);
    mem.B_sens = zeros(nx,nu*N);
    mem.Cx = zeros(nc,nx*N);
    mem.Cu = zeros(nc,nu*N);
    mem.gx = zeros(nx,N+1);
    mem.gu = zeros(nu,N);
    mem.a = zeros(nx,N);
    mem.ds0 = zeros(nx,1);
    mem.lc = zeros(N*nc+ncN,1);
    mem.uc = zeros(N*nc+ncN,1);
    mem.lb_du = zeros(N*nu,1);
    mem.ub_du = zeros(N*nu,1);
    mem.CxN = zeros(ncN,nx);

    mem.G = zeros(nx,N^2*nu);
    mem.Hc = zeros(N*nu,N*nu);
    mem.Cc = zeros(N*nc+ncN,N*nu);
    mem.gc = zeros(N*nu,1);
    mem.lcc = zeros(N*nc+ncN,1);
    mem.ucc = zeros(N*nc+ncN,1);

    mem.dz = zeros(nx+nu,N);
    mem.dxN= zeros(nx,1);
    mem.lambda_new = zeros(nx,N+1);
    mem.mu_new = zeros(nc,N);
    mem.muN_new = zeros(ncN,1);
    mem.mu_u_new = zeros(N*nu,1);
    
    if strcmp(opt.lin_obj,'yes')
        mem.lin_obj = 1;
        
        [Jx, Ju] = Ji_fun('Ji_fun',zeros(nx+nu,1),zeros(np,1),zeros(ny,1), input.W);
        Qi = full(Jx'*Jx);
        Si = full(Jx'*Ju);
        Ri = full(Ju'*Ju);
        mem.Q_h = repmat(Qi,1,N+1);
        mem.S = repmat(Si,1,N);
        mem.R = repmat(Ri,1,N);
        
        JN = JN_fun('JN_fun',zeros(nx,1),zeros(np,1),zeros(nyN,1), input.WN);
        mem.Q_h(:,N*nx+1:end) = full(JN'*JN);
    else
        mem.lin_obj = 0;
        mem.Q_h = zeros(nx,nx*(N+1));
        mem.S = zeros(nx,nu*N);
        mem.R = zeros(nu,nu*N);
    end
      
    %% input
    input.lambda=zeros(nx,N+1);
    input.mu=zeros(nc,N);
    input.muN=zeros(ncN,1);
    input.mu_u = zeros(N*nu,1);
    
end


function [mem] = InitMemory(settings, opt, input)

    Ts  = settings.Ts;       % Sampling time
    Ts_st = settings.Ts_st;  % Shooting interval
    s = settings.s;      % number of integration steps per interval
    nx = settings.nx;    % No. of states
    nu = settings.nu;    % No. of controls
    ny = settings.ny;    % No. of outputs (references)    
    nyN= settings.nyN;   % No. of outputs at terminal stage 
    np = settings.np;    % No. of parameters (on-line data)
    nbx = settings.nbx;  % No. of bounds on states
    nbu = settings.nbu;  % No. of bounds on controls
    nbx_idx = settings.nbx_idx; % indexes of states which are bounded
    nbu_idx = settings.nbu_idx; % indexes of controls which are bounded
    nc = settings.nc;    % No. of general constraints
    ncN = settings.ncN;  % No. of general constraints at terminal stage
    N     = settings.N;             % No. of shooting points
    N2  = settings.N2;

    %% memory
    mem = struct;
    mem.warm_start=0;
    mem.hot_start=0;
    if strcmp(opt.hotstart, 'yes')
        mem.hot_start=1;
    end
    
    if strcmp(opt.condensing,'partial_condensing')
        mem.settings2.N = N2;
        mem.settings2.nx = nx;
        mem.settings2.Nc = N/N2;
        mem.settings2.nu = mem.settings2.Nc*nu;
        mem.settings2.nbx = nbx;
        mem.settings2.nbx_idx = nbx_idx;
        
        mem.settings2.nbu_idx=[];
        for i=1:mem.settings2.Nc        
           mem.settings2.nbu_idx = [mem.settings2.nbu_idx,nbu_idx+(i-1)*nu];
        end
        
        mem.settings2.nc = nc*mem.settings2.Nc+(mem.settings2.Nc-1)*nbx;
        mem.settings2.ncN = ncN;
        
    end
    
    switch opt.qpsolver
        case 'qpoases'   
            mem.qpoases_opt = qpOASES_options('MPC');
%             mem.qpoases_opt = qpOASES_options('default');
        case 'qore'
            mem.qore_id = -1;
        case 'quadprog_dense'
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
        case 'ipopt_dense'
            ipopt_opts=ipoptset('constr_viol_tol',1e-3,'acceptable_tol',1e-3,'hessian_constant','yes',...
                        'mehrotra_algorithm','yes','mu_oracle','probing','jac_c_constant','yes',...
                        'jac_d_constant','yes','mu_strategy','adaptive','adaptive_mu_globalization',...
                        'never-monotone-mode','accept_every_trial_step','yes');
    
            mem.ipopt.options.eq=false(N*nc+ncN+N*nbx,1);
            mem.ipopt.options.ineq=true(N*nc+ncN+N*nbx,1);
            mem.ipopt.x0=zeros(N*nu,1);
            
            mem.ipopt.options.nleq=[];
            mem.ipopt.options.nlineq=[];
            mem.ipopt.options.ipopt=ipopt_opts;
            mem.ipopt.options.ipopt.print_level=0;
            
            
        case 'ipopt_sparse'
            nw = (N+1)*nx+N*nu;
            neq = (N+1)*nx;
            nineq = N*nc+ncN;
            
            ipopt_opts=ipoptset('constr_viol_tol',1e-3,'acceptable_tol',1e-3,...
                        'hessian_constant','yes','jac_c_constant','yes','jac_d_constant','yes',...
                        'mehrotra_algorithm','yes',...
                        'mu_oracle','probing',...
                        'mu_strategy','adaptive',...
                        'adaptive_mu_globalization','never-monotone-mode',...
                        'accept_every_trial_step','yes');
    
            mem.ipopt.options.eq=[false(nineq,1);true(neq,1)];
            mem.ipopt.options.ineq=[true(nineq,1);false(neq,1)];
            mem.ipopt.x0=zeros(nw,1);
            
            mem.ipopt.options.nleq=[];
            mem.ipopt.options.nlineq=[];
            mem.ipopt.options.ipopt=ipopt_opts;
            mem.ipopt.options.ipopt.print_level=0;
            
            mem.ipopt_data.H = zeros(nw,nw);
            mem.ipopt_data.g = zeros(nw,1);
            mem.ipopt_data.dG = zeros(neq,nw);  mem.ipopt_data.dG(1:nx,1:nx) = eye(nx);
            mem.ipopt_data.dBg = zeros(nineq,nw);
            mem.ipopt_data.G = zeros(neq,1);
            mem.ipopt.options.ub = inf(nw,1);
            mem.ipopt.options.lb = -inf(nw,1);
            
        case 'ipopt_partial_sparse'
                                
            nw = (N2+1)*mem.settings2.nx+N2*mem.settings2.nu;
            neq = (N2+1)*mem.settings2.nx;
            nineq = N2*mem.settings2.nc+mem.settings2.ncN;

            ipopt_opts=ipoptset('constr_viol_tol',1e-3,'acceptable_tol',1e-3,'hessian_constant','yes',...
                                'mehrotra_algorithm','yes','mu_oracle','probing','jac_c_constant','yes',...
                                'jac_d_constant','yes','mu_strategy','adaptive','adaptive_mu_globalization',...
                                'never-monotone-mode','accept_every_trial_step','yes');

            mem.mem2.ipopt.options.eq=[false(nineq,1);true(neq,1)];
            mem.mem2.ipopt.options.ineq=[true(nineq,1);false(neq,1)];
            mem.mem2.ipopt.x0=zeros(nw,1);

            mem.mem2.ipopt.options.nleq=[];
            mem.mem2.ipopt.options.nlineq=[];
            mem.mem2.ipopt.options.ipopt=ipopt_opts;
            mem.mem2.ipopt.options.ipopt.print_level=0;

            mem.mem2.ipopt_data.H = zeros(nw,nw);
            mem.mem2.ipopt_data.g = zeros(nw,1);
            mem.mem2.ipopt_data.dG = zeros(neq,nw);  mem.mem2.ipopt_data.dG(1:mem.settings2.nx,1:mem.settings2.nx) = eye(mem.settings2.nx);
            mem.mem2.ipopt_data.dBg = zeros(nineq,nw);
            mem.mem2.ipopt_data.G = zeros(neq,1);
            mem.mem2.ipopt.options.ub = inf(nw,1);
            mem.mem2.ipopt.options.lb = -inf(nw,1);

            
        case 'qpdunes'
            nw = (N+1)*nx+N*nu;
            mem.dunes.H=zeros(nx+nu,(nx+nu)*N);
            mem.dunes.P=zeros(nx,nx);
            mem.dunes.C=zeros(nx,(nx+nu)*N);
            mem.dunes.c=zeros(nx, N);
            mem.dunes.D=zeros(nc,(nx+nu)*N+nx);
            mem.dunes.g=zeros(nw,1);
            mem.dunes.zLow = -inf(nw,1);
            mem.dunes.zUpp = inf(nw,1);
            mem.dunes.dLow = -inf((N+1)*nc,1);
            mem.dunes.dUpp = inf((N+1)*nc,1);
            
            mem.dunes.qpOptions = qpDUNES_options( 'default', ...
                             'maxIter', 100, ...
                             'printLevel', 0, ...
                             'logLevel', 0, ...     % log all data
                             'lsType', 4, ...       % Accelerated gradient biscection LS
                             'stationarityTolerance', 1.e-6, ...
                             'regType', 2, ...       % regularize only singular directions; 1 is normalized Levenberg Marquardt
                             'newtonHessDiagRegTolerance', 1.e-8, ...
                             'maxNumLineSearchIterations',			19, ...			% 0.3^19 = 1e-10
                             'maxNumLineSearchRefinementIterations',	49, ...			% 0.62^49 = 1e-10
                             'lineSearchReductionFactor',		0.3, ...		% needs to be between 0 and 1
                             'lineSearchIncreaseFactor',			1.5 ...		% needs to be greater than 1
                             );
                         
        case 'osqp'
            nw = (N+1)*nx+N*nu;
            neq = (N+1)*nx;
            nineq = N*nc+ncN;
            
            mem.qp_obj = osqp;
            mem.osqp_options = mem.qp_obj.default_settings();
%             mem.osqp_options.eps_abs=1e-4;
%             mem.osqp_options.eps_rel=1e-4;
%             mem.osqp_options.polish = true;
            mem.osqp_options.verbose = false;
%             mem.osqp_options.linsys_solver = 'mkl pardiso';
            
            mem.osqp_data.H = zeros(nw,nw);
            mem.osqp_data.g = zeros(nw,1);
            mem.osqp_data.dG = zeros(neq,nw);  
            mem.osqp_data.dG(1:nx,1:nx) = eye(nx);
            mem.osqp_data.dBg = zeros(nineq,nw);
            mem.osqp_data.dBx = zeros(N*nbx,nw);
            mem.osqp_data.dBu = zeros(N*nu,nw);
            mem.osqp_data.G = zeros(neq,1);
            mem.osqp_data.ub = inf(nineq+N*nbx+N*nu,1);
            mem.osqp_data.lb = -inf(nineq+N*nbx+N*nu,1);
            
            for i=0:N-1
                for j=1:nbx
                    mem.osqp_data.dBx(i*nbx+1:(i+1)*nbx, (i+1)*nx+nbx_idx(j)) = 1;
                end
            end
            mem.osqp_data.dBu(:,neq+1:end) = eye(N*nu,N*nu);
                                     
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
            mem.newton_iter = 5;
            mem.JFK = mem.h*[mem.B_tab(1)*eye(nx,nx), mem.B_tab(2)*eye(nx,nx), mem.B_tab(3)*eye(nx,nx)];
        otherwise 
            error('Please choose a correct integrator');       
    end
    
    % globalization
    mem.sqp_maxit = 1;           % maximum number of iterations for each sampling instant (for RTI, this is ONE)
    mem.kkt_lim = 1e-4;          % tolerance on optimality
    mem.mu_merit=0;              % initialize the parameter
    mem.eta=1e-4;                % merit function parameter
    mem.tau=0.8;                 % step length damping factor
    mem.mu_safty=1.1;            % constraint weight update factor (for merit function)
    mem.rho=0.5;                 % merit function parameter
    mem.alpha=1;                 % default step length
    
    % allocate memory
    mem.A = zeros(nx,nx*N);
    mem.B = zeros(nx,nu*N);
    mem.Cx = zeros(nbx,nx);
    mem.Cgx = zeros(nc,nx*N);
    mem.Cgu = zeros(nc,nu*N);
    mem.CgN = zeros(ncN,nx);
    mem.gx = zeros(nx,N+1);
    mem.gu = zeros(nu,N);
    mem.a = zeros(nx,N);
    mem.ds0 = zeros(nx,1);
    mem.lc = zeros(N*nc+ncN,1);
    mem.uc = zeros(N*nc+ncN,1);
    mem.lb_du = zeros(N*nu,1);
    mem.ub_du = zeros(N*nu,1);
    mem.lb_dx = zeros(N*nbx,1);
    mem.ub_dx = zeros(N*nbx,1);
    
    mem.Hc = zeros(N*nu,N*nu);
    mem.Ccx = zeros(N*nbx,N*nu);
    mem.Ccg = zeros(N*nc+ncN,N*nu);
    mem.gc = zeros(N*nu,1);
    mem.lcc = zeros(N*nc+ncN,1);
    mem.ucc = zeros(N*nc+ncN,1);
    mem.lxc = zeros(N*nbx,1);
    mem.uxc = zeros(N*nbx,1);
    
    mem.dx = zeros(nx,N+1);
    mem.du = zeros(nu,N);
    mem.lambda_new = zeros(nx,N+1);
    mem.mu_new = zeros(N*nc+ncN,1);
    mem.mu_x_new = zeros(N*nbx,1);
    mem.mu_u_new = zeros(N*nu,1);
    
    mem.q_dual = zeros(nx,N+1);
    mem.dmu = zeros(N*nu+N*nbx+N*nc+ncN,1);
    
    for i=1:nbx
        mem.Cx(i,nbx_idx(i)) = 1.0;
    end
    
    if strcmp(opt.lin_obj,'yes')
        mem.lin_obj = 1;
        
        [Jx, Ju] = Ji_fun('Ji_fun',zeros(nx,1),zeros(nu,1),zeros(np,1),zeros(ny,1), input.W);
        Qi = full(Jx'*Jx);
        for i=1:nx
            if Qi(i,i)<1e-12
                Qi(i,i)=1e-12;
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
            if Qf(i,i)<1e-12
                Qf(i,i)=1e-12;
            end
        end
        mem.Q(:,N*nx+1:end) = Qf;
    else
        mem.lin_obj = 0;
        mem.Q = zeros(nx,nx*(N+1));
        mem.S = zeros(nx,nu*N);
        mem.R = zeros(nu,nu*N);
    end
    mem.reg = 1e-12;
              
    mem.iter=1;
    
     %% for CMON-RTI	
    mem.F_old = zeros(nx,N);	
    mem.CMON_pri = zeros(N,1);	
    mem.CMON_dual = zeros(N,1);	
    mem.q_dual = zeros(nx,N+1);	
    mem.V_pri = zeros(nx,N);
    mem.V_dual = zeros(nx+nu,N);
    mem.dmu = zeros(N*nu+N*nbx+N*nc+ncN,1);
    mem.threshold_pri = 0;	
    mem.threshold_dual = 0;	
    mem.tol=0;	
    mem.perc=100;
    mem.idxc=zeros(N,1);
    
    mem.tol_abs=5e-1;
    mem.tol_ref=5e-1;  	       
    mem.alpha_cmon = 1;      
    mem.beta_cmon = 1;        
    mem.c1 = 0.1;
    mem.gamma = 0;	
    mem.rho_cmon = 0;
          
    mem.local = 0;
    
    mem.rho_ratio=[];
    mem.gamma_ratio=[];
    
    mem.r_ratio=[];
    
    mem.shift_x = zeros(nx,N+1);
    mem.shift_u = zeros(nu,N);
    
    mem.Ns=25;
       
end


clear mex; close all;clc;
%% Configuration (complete your configuration here...)
addpath('nmpc');
addpath('Source_Codes');
addpath('mex_core');

if exist('settings','file')==2
    load settings
else 
    error('No setting data is detected!');
end

Ts  = settings.Ts;   % Sampling time
nx = settings.nx;    % No. of states
nu = settings.nu;    % No. of controls
ny = settings.ny;    % No. of outputs (references)    
nyN= settings.nyN;   % No. of outputs at terminal stage 
np = settings.np;    % No. of parameters (on-line data)
nc = settings.nc;    % No. of constraints
ncN = settings.ncN;  % No. of constraints at terminal stage

N     = 20; % No. of shooting points
nw    = (N+1)*nx+N*nu; % No. of total optimization varialbes
neq   = (N+1)*nx; % No. of equality constraints
nineq = N*2*nc+2*ncN; % No. of inequality constraints (by default we assume there are lower and upper bounds)

settings.N  = N;
settings.nw = nw; 
settings.neq = neq;
settings.nineq = nineq; 

% solver configurations
input.opt.integrator='ERK4'; % 'ERK4', 'IRK4'
input.opt.hessian='gauss_newton';  % 'gauss_newton', 'exact'
input.opt.qpsolver='qpoases'; %'qpoases'
input.opt.condensing='full';  %'full','no','null-space'
input.opt.hotstart='no'; %'yes','no' (only for qpoases)
input.opt.shifting='no'; % 'yes','no'

% globalization configurations (require advanced knowledge on globilization)
% input.opt.meritfun.mu_merit=0;              % initialize the parameter
% input.opt.meritfun.eta=1e-4;                % merit function parameter
% input.opt.meritfun.tau=0.8;                 % step length damping factor
% input.opt.meritfun.mu_safty=1.1;            % constraint weight update factor (for merit function)
% input.opt.meritfun.rho=0.5;                 % merit function parameter
% input.opt.meritfun.second_order_correction='no';  % trigger on second_order_correction

%% Initialize all solvers (skip if you don't understand, we will take care of everything)

% if strcmp(input.opt.qpsolver,'ipopt')
%     ipopt_opts=ipoptset('constr_viol_tol',1e-3,'acceptable_tol',1e-3,'hessian_constant','yes',...
%                         'mehrotra_algorithm','yes','mu_oracle','probing','jac_c_constant','yes',...
%                         'jac_d_constant','yes','mu_strategy','adaptive','adaptive_mu_globalization',...
%                         'never-monotone-mode','accept_every_trial_step','yes');
%     if strcmp(input.opt.condensing,'no')
%         input.opt.ipopt.options.eq=[false(nineq,1);true(neq,1)];
%         input.opt.ipopt.options.ineq=[true(nineq,1);false(neq,1)];
%         input.opt.ipopt.x0=zeros(nw,1);
%     else
%         input.opt.ipopt.options.eq=false(settings.nineq,1);
%         input.opt.ipopt.options.ineq=true(settings.nineq,1);
%         input.opt.ipopt.x0=zeros(N*nu,1);
%     end
%     input.opt.ipopt.options.nleq=[];
%     input.opt.ipopt.options.nlineq=[];
%     input.opt.ipopt.options.ipopt=ipopt_opts;
%     input.opt.ipopt.options.ipopt.print_level=0;
% end

if strcmp(input.opt.qpsolver,'qpoases')
    input.opt.condensing_matrix.QP=0;
end


%% Off-line data structure initialization (skip if you don't understand)

% Multipliers

input.lambda=ones(nx,N+1);
input.mu=zeros(2*nc,N);
input.muN=zeros(2*ncN,1);

% ERK settings

sim_erk_mem.A_erk=[0, 0, 0, 0;
               0.5, 0, 0, 0;
               0, 0.5, 0, 0;
               0, 0, 1, 0];
sim_erk_mem.B_erk=[1/6, 1/3, 1/3, 1/6];
sim_erk_mem.num_steps = 2;
sim_erk_mem.num_stages = 4;
sim_erk_mem.h=0.1/sim_erk_mem.num_steps;
sim_erk_mem.nx = nx;
sim_erk_mem.nu = nu;
sim_erk_mem.Sx = eye(nx);
sim_erk_mem.Su = zeros(nx,nu);

%% IRK settings
sim_irk_mem.A_irk=[5/36,             2/9-sqrt(15)/15, 5/36-sqrt(15)/30;
                   5/36+sqrt(15)/24, 2/9            , 5/36-sqrt(15)/24;
                   5/36+sqrt(15)/30, 2/9+sqrt(15)/15, 5/36];
sim_irk_mem.B_irk=[5/18;4/9;5/18];
sim_irk_mem.num_steps = 2;
sim_irk_mem.num_stages = 3;
sim_irk_mem.h= 0.1/sim_irk_mem.num_steps;
sim_irk_mem.nx = nx;
sim_irk_mem.nu = nu;
sim_irk_mem.Sx = eye(nx);
sim_irk_mem.Su = zeros(nx,nu);
sim_irk_mem.newton_iter = 3;
sim_irk_mem.JFK = sim_irk_mem.h*[sim_irk_mem.B_irk(1)*eye(nx,nx), sim_irk_mem.B_irk(2)*eye(nx,nx), sim_irk_mem.B_irk(3)*eye(nx,nx)];
%% Initialzation (initialize your simulation properly...)

Initialization;

% memory_init(sim_irk_mem);
%% Simulation (start your simulation...)

iter = 1; time = 0.0;
Tf = 4;               % simulation time
state_sim= x0';
controls_MPC = u0';
y_sim = [];
constraints = [];
CPT = [];

% mpc_callnum=1;         % maximum number of iterations for each sampling instant (for RTI, this is ONE)
% kkt_lim=1e-2;           % tolerance on optimality

while time(end) < Tf
    
    % the reference input.y is a ny by N matrix
    % the reference input.yN is a nyN by 1 vector
    
    % time-invariant reference
    input.y = repmat(REF',1,N);
    input.yN = REF(1:nyN)';
    
    % time-varying reference (no reference preview)
%     input.y = repmat(REF(iter,:)',1,N);
%     input.yN = REF(iter,1:nyN)';

    % time-varying reference (reference preview)
%     input.y = REF(iter:iter+N-1,:);
%     input.yN = REF(iter+N,:);
           
    % obtain the state measurement
    input.x0 = state_sim(end,:)';
    
    % call the NMPC solver
    output=mpc_nmpcsolver(input,settings,sim_erk_mem,sim_irk_mem);
    
    % obtain the solution and update the data
    switch input.opt.shifting
        case 'yes'
        input.z=[output.z(:,2:end),[output.xN; output.z(nx+1:nx+nu,end)]];  
        input.xN=output.xN;
        input.lambda=[output.lambda(:,2:end),output.lambda(:,end)];
        input.mu=[output.mu(:,2:end),[output.muN;output.mu(2*ncN+1:2*nc,end)]];
        input.muN=output.muN;
        case 'no'
        input.z=output.z;
        input.xN=output.xN;
        input.lambda=output.lambda;
        input.mu=output.mu;
        input.muN=output.muN;
    end
%     input.opt.meritfun=output.meritfun;
    
    if strcmp(input.opt.qpsolver,'qpoases')
        input.opt.condensing_matrix.QP=output.opt.condensing_matrix.QP;
    end
    
    % collect the statistics
    cpt=output.info.cpuTime;
    tshooting=output.info.shootTime;
    tcond=output.info.condTime;
    tqp=output.info.qpTime;
    KKT=output.info.kktValue;
    
    % Simulate system dynamics
    sim_input.x = state_sim(end,:).';
    sim_input.u = output.z(nx+1:nx+nu,1);
    sim_input.p = input.od(:,1)';
    xf=full( Simulate_system('Simulate_system', sim_input.x, sim_input.u, sim_input.p) ); 
    
    % Collect outputs
    y_sim = [y_sim; full(h_fun('h_fun', xf, sim_input.u, sim_input.p))'];  
    
    % Collect constraints
    constraints=[constraints; full( ineq_fun('ineq_fun', xf, sim_input.u, sim_input.p) )'];
    
    % store the optimal solution and states
    controls_MPC = [controls_MPC; output.z(nx+1:nx+nu,1)'];
    state_sim = [state_sim; xf'];
    
    % go to the next sampling instant
    nextTime = iter*Ts; 
    iter = iter+1;
    disp(['current time:' num2str(nextTime) '  CPT:' num2str(cpt) 'ms  MULTIPLE SHOOTING:' num2str(tshooting) 'ms  COND:' num2str(tcond) 'ms  QP:' num2str(tqp) 'ms  KKT:' num2str(KKT)]);
    time = [time nextTime];
    
    CPT = [CPT; cpt, tshooting, tcond, tqp];
end

%% draw pictures (optional)
clear mex;

mean(CPT,1)
Draw;

clear mex; close all;clc;

%% Configuration (complete your configuration here...)
addpath([pwd,'/nmpc']);
addpath([pwd,'/Source_Codes']);
addpath([pwd,'/mex_core']);

if exist('settings','file')==2
    load settings
else 
    error('No setting data is detected!');
end

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

%% solver configurations
N  = 40;             % No. of shooting points
settings.N = N;

opt.integrator='ERK4-CASADI'; % 'ERK4','IRK3, 'ERK4-CASADI'(for test)
opt.hessian='gauss_newton';  % 'gauss_newton', 'exact'
opt.qpsolver='qpoases'; %'qpoases'
opt.condensing='full';  %'full'
opt.hotstart='no'; %'yes','no' (only for qpoases)
opt.shifting='no'; % 'yes','no'

%% Initialize Solvers (only for advanced users)

[input, mem] = InitMemory(settings, opt, input);

[input_exact, mem_exact] = InitMemory(settings, opt, input_exact);

%% Initialize Data (all users have to do this)

[input, data] = InitData(settings, input);

[input_exact, data_exact] = InitData(settings, input_exact);

%% Simulation (start your simulation...)

iter = 1; time = 0.0;
Tf = 50;               % simulation time
state_sim= x0';
controls_MPC = u0';
y_sim = [];
constraints = [];
CPT = [];
ERR = [];

ref_traj = [];

while time(end) < Tf
    
    % the reference input.y is a ny by N matrix
    % the reference input.yN is a nyN by 1 vector
    
    % time-invariant reference
    input.y = repmat(REF',1,N);
    input.yN = REF(1:nyN)';
    
    % time-varying reference (no reference preview)
%     input.y = repmat(REF(iter,:)',1,N);
%     input.yN = REF(iter,1:nyN)';
    
    %time-varying reference (reference preview)
%     REF = zeros(ny,N+1);
%     for i=1:N+1
%         x = amplitude_x*sin(((time(end)+(i-1)*Ts_st))*2*pi*f_x);
%         theta = amplitude_theta*sin(((time(end)+(i-1)*Ts_st))*2*pi*f_theta);
%         REF(:,i) = [x 0 0 0 theta 0 zeros(1,nu)]';
%     end
%     ref_traj=[ref_traj, REF(:,1)];
%     input.y = REF(:,1:N);
%     input.yN = REF(1:nyN,N+1);
           
    % obtain the state measurement
    input.x0 = state_sim(end,:)';
    
    % call the NMPC solver
    input_exact.z(:) = input.z(:);
    input_exact.xN(:) = input.xN(:);
    input_exact.lambda(:) = input.lambda(:);
    input_exact.y(:) = input.y(:);
    input_exact.yN = input.yN(:);
    input_exact.x0 = input.x0(:);
    
    [output, mem]=mpc_nmpcsolver(input, settings, mem, 0);
    dw = [reshape(mem.dz, settings.N*(settings.nx+settings.nu), 1); mem.dxN];
    dl = reshape(mem.q_dual, (settings.N+1)*settings.nx, 1);
    dm = input.mu_u;
    dy = [dw;dm;dl];
    
    [output_exact, mem_exact]=mpc_nmpcsolver(input_exact, settings, mem_exact, 1);
    dw_exact = [reshape(mem_exact.dz, settings.N*(settings.nx+settings.nu), 1); mem_exact.dxN];
    dl_exact = reshape(mem_exact.q_dual, (settings.N+1)*settings.nx, 1);
    dm_exact = input_exact.mu_u;
    dy_exact = [dw_exact;dm_exact;dl_exact];
        
    ERR = [ERR;norm(dy-dy_exact)];
    
    % obtain the solution and update the data
    switch opt.shifting
        case 'yes'
        input.z=[output.z(:,2:end),[output.xN; output.z(nx+1:nx+nu,end)]];  
        input.xN=output.xN;
        input.lambda=[output.lambda(:,2:end),output.lambda(:,end)];
        input.mu=[output.mu(:,2:end),[output.muN;output.mu(ncN+1:nc,end)]];
        input.muN=output.muN;
        case 'no'
        input.z=output.z;
        input.xN=output.xN;
        input.lambda=output.lambda;
        input.mu=output.mu;
        input.muN=output.muN;
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
    constraints=[constraints; full( path_con_fun('path_con_fun', xf, sim_input.u, sim_input.p) )'];
    
    % store the optimal solution and states
    controls_MPC = [controls_MPC; output.z(nx+1:nx+nu,1)'];
    state_sim = [state_sim; xf'];
    
    % go to the next sampling instant
    nextTime = iter*Ts; 
    iter = iter+1;
    disp(['current time:' num2str(nextTime) '  CPT:' num2str(cpt) 'ms  MULTIPLE SHOOTING:' num2str(tshooting) 'ms  COND:' num2str(tcond) 'ms  QP:' num2str(tqp) 'ms  KKT:' num2str(KKT)]);
    disp(['exactly updated sensitivities:' num2str(mem.perc) '%']);
    disp(['threshold_pri:' num2str(mem.threshold_pri) '   threshold_dual:' num2str(mem.threshold_dual) '   ERR:' num2str(ERR(end)) '   Tol:' num2str(mem.tol)]);
%     disp(['No. of SQP Iteration: ' num2str(output.info.iteration_num)]);
    disp('   ');
    
    time = [time nextTime];
    
    CPT = [CPT; cpt, tshooting, tcond, tqp];
end

qpOASES_sequence( 'c', mem.qpoases.warm_start);
clear mex;

%% draw pictures (optional)
disp('Average CPT:');
mean(CPT(2:end-1,:),1)

disp('Maximum CPT: ')
max(CPT(2:end-1,:))

Draw;

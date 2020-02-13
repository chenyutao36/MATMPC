clear all; clear mex; close all;clc;

disp( ' ' );
disp( 'MATMPC -- A (MAT)LAB based Model(M) Predictive(P) Control(C) Package.' );
disp( 'Copyright (C) 2016-2019 by Yutao Chen, University of Padova' );
disp( 'All rights reserved.' );
disp( ' ' );
disp( 'MATMPC is distributed under the terms of the' );
disp( 'GNU General Public License 3.0 in the hope that it will be' );
disp( 'useful, but WITHOUT ANY WARRANTY; without even the implied warranty' );
disp( 'of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.' );
disp( 'See the GNU General Public License for more details.' );
disp( ' ' );
disp( ' ' );
disp('---------------------------------------------------------------------------------');

%% Configuration (complete your configuration here...)
addpath([pwd,'/nmpc']);
addpath([pwd,'/model_src']);
addpath([pwd,'/mex_core']);
addpath(genpath([pwd,'/data']));
if ismac
    addpath(genpath([pwd,'/solver/mac']));
elseif isunix
    addpath(genpath([pwd,'/solver/linux']));
elseif ispc
    addpath(genpath([pwd,'/solver/win64']));
else
    disp('Platform not supported')
end

cd data;
if exist('settings','file')==2
    load settings
    cd ..
else 
    cd ..
    error('No setting data is detected!');
end

Ts = settings.Ts_st;     % Closed-loop sampling time (usually = shooting interval)

Ts_st = settings.Ts_st;  % Shooting interval
nx = settings.nx;    % No. of states
nu = settings.nu;    % No. of controls
ny = settings.ny;    % No. of outputs (references)    
nyN= settings.nyN;   % No. of outputs at terminal stage 
np = settings.np;    % No. of parameters (on-line data)
nc = settings.nc;    % No. of constraints
ncN = settings.ncN;  % No. of constraints at terminal stage
nbx = settings.nbx;  % No. of state bounds

%% solver configurations

N  = 80;             % No. of shooting points
settings.N = N;

N2 = N/5;
settings.N2 = N2;    % No. of horizon length after partial condensing (N2=1 means full condensing)

r = 10;
settings.r = r;      % No. of input blocks (go to InitMemory.m, line 418 to configure)

opt.hessian         = 'Gauss_Newton';  % 'Gauss_Newton', 'Generalized_Gauss_Newton'
opt.integrator      = 'ERK4'; % 'ERK4','IRK3','IRK3-DAE'
opt.condensing      = 'default_full';  %'default_full','no','blasfeo_full(require blasfeo installed)','partial_condensing'
opt.qpsolver        = 'qpoases'; 
opt.hotstart        = 'no'; %'yes','no' (only for qpoases, use 'no' for nonlinear systems)
opt.shifting        = 'no'; % 'yes','no'
opt.ref_type        = 0; % 0-time invariant, 1-time varying(no preview), 2-time varying (preview)
opt.nonuniform_grid = 0; % if use non-uniform grid discretization (go to InitMemory.m, line 436 to configure)
opt.RTI             = 'yes'; % if use Real-time Iteration
%% available qpsolver

%'qpoases' (condensing is needed)
%'qpoases_mb' (move blocking strategy)
%'quadprog_dense' (for full condensing)
%'hpipm_sparse' (run mex_core/compile_hpipm.m first; set opt.condensing='no')
%'hpipm_pcond' (run mex_core/compile_hpipm.m first; set opt.condensing='no')
%'ipopt_dense' (install OPTI Toolbox first; for full condensing)
%'ipopt_sparse' (install OPTI Toolbox first; set opt.condensing='no')
%'ipopt_partial_sparse'(set opt.condensing='partial_condensing'; only for state and control bounded problems)
%'osqp_sparse' (set opt.condensing='no')
%'osqp_partial_sparse' (set opt.condensing='partial_condensing')
%'qpalm_cond' (condensing is needed)
%'qpalm_sparse'(set opt.condensing='no')

%% Initialize Data (all users have to do this)
if opt.nonuniform_grid
    [input, data] = InitData_ngrid(settings);
    N = r;
    settings.N = N;
else
	[input, data] = InitData(settings);
end  

%% Initialize Solvers (only for advanced users)

mem = InitMemory(settings, opt, input);

%% Simulation (start your simulation...)

mem.iter = 1; time = 0.0;
Tf = 4;  % simulation time
state_sim= input.x0';
controls_MPC = input.u0';
y_sim = [];
constraints = [];
CPT = [];
ref_traj = [];
KKT = [];
OBJ=[];
numIT=[];

while time(end) < Tf
        
    % the reference input.y is a ny by N matrix
    % the reference input.yN is a nyN by 1 vector    
    switch opt.ref_type
        case 0 % time-invariant reference
            input.y = repmat(data.REF',1,N);
            input.yN = data.REF(1:nyN)';
        case 1 % time-varying reference (no reference preview)
            input.y = repmat(data.REF(mem.iter,:)',1,N);
            input.yN = data.REF(mem.iter,1:nyN)';
        case 2 %time-varying reference (reference preview)
            input.y = data.REF(mem.iter:mem.iter+N-1,:)';
            input.yN = data.REF(mem.iter+N,1:nyN)';
    end
              
    % obtain the state measurement
    input.x0 = state_sim(end,:)';
    
    % call the NMPC solver 
    [output, mem] = mpc_nmpcsolver(input, settings, mem, opt);
        
    % obtain the solution and update the data
    switch opt.shifting
        case 'yes'
            input.x=[output.x(:,2:end),output.x(:,end)];  
            input.u=[output.u(:,2:end),output.u(:,end)];
            input.z=[output.z(:,2:end),output.z(:,end)];
            input.lambda=[output.lambda(:,2:end),output.lambda(:,end)];
            input.mu=[output.mu(nc+1:end);output.mu(end-nc+1:end)];
            input.mu_x=[output.mu_x(nbx+1:end);output.mu_x(end-nbx+1:end)];
            input.mu_u=[output.mu_u(nu+1:end);output.mu_u(end-nu+1:end)];
            
            % for CMoN-RTI
%             mem.A=[mem.A(:,nx+1:end),mem.A(:,end-nx+1:end)];
%             mem.B=[mem.B(:,nu+1:end),mem.B(:,end-nu+1:end)];
%             mem.F_old = [mem.F_old(:,2:end),mem.F_old(:,end)];
%             mem.V_pri = [mem.V_pri(:,2:end),mem.V_pri(:,end)];
%             mem.V_dual = [mem.V_dual(:,2:end),mem.V_dual(:,end)];
%             mem.q_dual = [mem.q_dual(:,2:end),mem.q_dual(:,end)];            
%             mem.shift_x = input.x-output.x;
%             mem.shift_u = input.u-output.u;
        case 'no'
            input.x=output.x;
            input.u=output.u;
            input.z=output.z;
            input.lambda=output.lambda;
            input.mu=output.mu;
            input.mu_x=output.mu_x;
            input.mu_u=output.mu_u;
    end
    
    % collect the statistics
    cpt=output.info.cpuTime;
    tshooting=output.info.shootTime;
    tcond=output.info.condTime;
    tqp=output.info.qpTime;
    OptCrit=output.info.OptCrit;
    
    % Simulate system dynamics
    sim_input.x = state_sim(end,:).';
    sim_input.u = output.u(:,1);
    sim_input.z = input.z(:,1);
    sim_input.p = input.od(:,1);

    [xf, zf] = Simulate_System(sim_input.x, sim_input.u, sim_input.z, sim_input.p, mem, settings);
    xf = full(xf);
    
    % Collect outputs
    y_sim = [y_sim; full(h_fun('h_fun', xf, sim_input.u, sim_input.p))'];  
    
    % Collect constraints
    constraints=[constraints; full( path_con_fun('path_con_fun', xf, sim_input.u, sim_input.p) )'];
        
    % store the optimal solution and states
    controls_MPC = [controls_MPC; output.u(:,1)'];
    state_sim = [state_sim; xf'];
    KKT= [KKT;OptCrit];
    OBJ= [OBJ;output.info.objValue];
    CPT = [CPT; cpt, tshooting, tcond, tqp];
    numIT = [numIT; output.info.iteration_num];
    
    % go to the next sampling instant
    nextTime = mem.iter*Ts; 
    mem.iter = mem.iter+1;
    disp(['current time:' num2str(nextTime) '  CPT:' num2str(cpt) 'ms  SHOOTING:' num2str(tshooting) 'ms  COND:' num2str(tcond) 'ms  QP:' num2str(tqp) 'ms  Opt:' num2str(OptCrit) '   OBJ:' num2str(OBJ(end)) '  SQP_IT:' num2str(output.info.iteration_num)]);
%     disp(['current time:' num2str(nextTime) '  CPT:' num2str(cpt) 'ms  SHOOTING:' num2str(tshooting) 'ms  COND:' num2str(tcond) 'ms  QP:' num2str(tqp) 'ms  Opt:' num2str(OptCrit) '   OBJ:' num2str(OBJ(end)) '  SQP_IT:' num2str(output.info.iteration_num) '  Perc:' num2str(mem.perc)]);   
    time = [time nextTime];   
end

%%
if strcmp(opt.qpsolver, 'qpoases')
    qpOASES_sequence( 'c', mem.warm_start);
end
% if strcmp(opt.qpsolver, 'qpalm')
%     mem.qpalm_solver.delete();
% end
clear mex;

%% draw pictures (optional)
disp(['Average CPT: ', num2str(mean(CPT(2:end,:),1)) ]);
disp(['Maximum CPT: ', num2str(max(CPT(2:end,:))) ]);

Draw;

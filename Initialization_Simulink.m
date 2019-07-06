% ------------------------------------------------------
%  This is an example of initializing simulink simulation
%  ------------------------------------------------------

%%
clear mex; close all; clear; clc;

addpath([pwd,'/nmpc']);
addpath([pwd,'/model_src']);
addpath([pwd,'/mex_core']);
%% Parametri Simulazione
cd data;
if exist('settings','file')==2
    load('settings');
    cd ..
else 
    cd ..
    error('No setting data is detected!');
end

Ts  = settings.Ts;       % Sampling time
Ts_st = settings.Ts_st;  % Shooting interval
s = settings.s;          % number of integration steps per interval
nx = settings.nx;        % No. of states
nu = settings.nu;        % No. of controls
ny = settings.ny;        % No. of outputs (references)    
nyN= settings.nyN;       % No. of outputs at terminal stage 
np = settings.np;        % No. of parameters (on-line data)
nc = settings.nc;        % No. of constraints
ncN = settings.ncN;      % No. of constraints at terminal stage
nbu = settings.nbu;      % No. of control bounds
nbx = settings.nbx;      % No. of state bounds
nbu_idx = settings.nbu_idx; % Index of control bounds
nbx_idx = settings.nbx_idx; % Index of state bounds

%% add more to Settings

N  = 80;
N2 = 5;
r  = 10;

settings.N = N;
settings.N2 = N2;
settings.r = r;

%% options
opt.hessian='Gauss_Newton';  % 'Gauss_Newton', 'Generalized_Gauss_Newton'
opt.integrator='ERK4'; % 'ERK4','IRK3, 'ERK4-CASADI'
opt.condensing='default_full';  %'default_full','no','blasfeo_full(require blasfeo installed)','partial_condensing'
opt.qpsolver='qpoases'; 
opt.hotstart='no'; %'yes','no' (only for qpoases)
opt.shifting='no'; % 'yes','no'
opt.lin_obj='no'; % 'yes','no' % if the inner objective function is linear and the outer objective is sum of quadratic
opt.ref_type=0; % 0-time invariant, 1-time varying(no preview), 2-time varying (preview)
opt.nonuniform_grid=0; % currently not supported 

%% available qpsolver
%'qpoases' (for full condensing)
%'qpoases_mb' (for full condensing+moving block, please use ERK4 as the integrator)
%'hpipm_sparse' (run mex_core/compile_hpipm.m first; set opt.condensing='no')
%'hpipm_pcond' (run mex_core/compile_hpipm.m first; set opt.condensing='no')
 
%% Initialization

x0 = [0;pi;0;0];    
u0 = zeros(nu,1);    
para0 = zeros(max(1,np),1);  

W=repmat([10 10 0.1 0.1 0.01]',1,N);
WN=W(1:nyN,1);

% upper and lower bounds for states (=nbx)
lb_x = -2;
ub_x = 2;

% upper and lower bounds for controls (=nbu)           
lb_u = -20;
ub_u = 20;
                       
% upper and lower bounds for general constraints (=nc)
lb_g = [];
ub_g = [];            
lb_gN = [];
ub_gN = [];
 
lb = repmat(lb_g,N,1);
ub = repmat(ub_g,N,1);
lb = [lb;lb_gN];
ub = [ub;ub_gN];
if isempty(lb)
    lb=0;
    ub=0;
end
        
lbu = -inf(nu,1);
ubu = inf(nu,1);
for i=1:nbu
    lbu(nbu_idx(i)) = lb_u(i);
    ubu(nbu_idx(i)) = ub_u(i);
end
        
lbu = repmat(lbu,1,N);
ubu = repmat(ubu,1,N);

lbx = repmat(lb_x,1,N+1);
ubx = repmat(ub_x,1,N+1);
if isempty(lbx)
    lbx=0;
    ubx=0;
end

x = repmat(x0,1,N+1);   
u = repmat(u0,1,N);    
para = repmat(para0,1,N+1);  

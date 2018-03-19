clear mex; close all; clear; clc;
%% Configuration (complete your configuration here...)
addpath('/home/chen/Documents/Packages/MATMPC/nmpc');
addpath('/home/chen/Documents/Packages/MATMPC/model_src');
addpath('/home/chen/Documents/Packages/MATMPC/mex_core');
%% Parametri Simulazione

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
nbu = settings.nbu;
nbx = settings.nbx;
nbu_idx = settings.nbu_idx;
nbx_idx = settings.nbx_idx;

%% Prediction Horizon
N=40;
settings.N = N;

%%
% global x0 lb ub lbN ubN lbu ubu W WN para ref;

x0 = [0;pi;0;0];    
u0 = zeros(nu,1);    
para0 = 0;  

W=diag([10 10 0.1 0.1 0.01]);
WN=W(1:nyN,1:nyN);

% upper and lower bounds for states (=nbx)
lb_x = -2;
ub_x = 2;
lb_xN = -2;
ub_xN = 2;

% upper and lower bounds for controls (=nbu)           
lb_u = -20;
ub_u = 20;
                       
% upper and lower bounds for general constraints (=nc)
lb_g = [];
ub_g = [];            
lb_gN = [];
ub_gN = [];

lb=repmat([lb_g;lb_x],1,N);
ub=repmat([ub_g;ub_x],1,N); 
lbN=[lb_gN;lb_x];               
ubN=[ub_gN;ub_x]; 
        
lbu = -inf(nu,1);
ubu = inf(nu,1);
for i=1:nbu
    lbu(nbu_idx(i)) = lb_u(i);
    ubu(nbu_idx(i)) = ub_u(i);
end
        
lbu = repmat(lbu,1,N);
ubu = repmat(ubu,1,N);
x = repmat(x0,1,N+1);   
u = repmat(u0,1,N);    
para = repmat(para0,1,N+1);  

ref = zeros(ny,1);

%%
assignin('base','Ts',Ts);
assignin('base','x0',x0);
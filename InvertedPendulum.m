%------------------------------------------%
% Inverted Pendulum 
  
% from "Autogenerating microsecond solvers for nonlinear MPC: A tutorial
% using ACADO integrators", Quirynen, 2015

%------------------------------------------%


%% Dimensions

nx=4;  % No. of states
nu=1;  % No. of controls
ny=5; % No. of outputs
nyN=4; % No. of outputs at the terminal point
np=0; % No. of model parameters
nc=0; % No. of general constraints
ncN=0; % No. of general constraints at the terminal point
nbx = 1; % No. of bounds on states
nbu = 1; % No. of bounds on controls


import casadi.*

states   = SX.sym('states',nx,1);
controls = SX.sym('controls',nu,1);
params   = SX.sym('paras',np,1);
refs     = SX.sym('refs',ny,1);
refN     = SX.sym('refs',nyN,1);
Q        = SX.sym('Q',ny,ny);
QN       = SX.sym('QN',nyN,nyN);


%% Dynamics

M = 1; 
m = 0.1;
l = 0.8; 
g = 9.81;

p=states(1);
theta=states(2);
v=states(3);
omega=states(4);  
u=controls(1);

a=-m*l*sin(theta)*omega^2+m*g*cos(theta)*sin(theta)+u;
b=-m*l*cos(theta)*omega^2+u*cos(theta)+(M+m)*g*sin(theta);
c=M+m-m*(cos(theta))^2;

x_dot=[v;omega;a/c;b/(l*c)];
xdot = SX.sym('xdot',nx,1);
impl_f = xdot - x_dot;
     
%% Objectives and constraints

% objectives
h = [p;theta;v;omega;u];

hN = h(1:nyN);

h_fun=Function('h_fun', {states,controls,params}, {h},{'states','controls','params'},{'h'});
hN_fun=Function('hN_fun', {states,params}, {hN},{'states','params'},{'hN'});

% general inequality constraints
general_con = [];
general_con_N = [];

% state and control bounds
nbx_idx = 1; % indexs of states which are bounded
nbu_idx = 1; % indexs of controls which are bounded
path_con=general_con;
path_con_N=general_con_N;
for i=1:nbx
    path_con=[path_con;states(nbx_idx(i))];
    path_con_N=[path_con_N;states(nbx_idx(i))];
end    
nc=nc+nbx;
ncN=ncN+nbx;

% build the function for inequality constraints
path_con_fun=Function('path_con_fun', {states,controls,params}, {path_con},{'states','controls','params'},{'path_con'});
path_con_N_fun=Function('path_con_N_fun', {states,params}, {path_con_N},{'states','params'},{'path_con_N'});


%% NMPC sampling time [s]

Ts = 0.05; % simulation sample time
Ts_st = 0.05; % shooting interval time
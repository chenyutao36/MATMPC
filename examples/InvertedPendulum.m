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

% state and control bounds
nbx_idx = 1; % indexs of states which are bounded
nbu_idx = 1; % indexs of controls which are bounded

%% create variables

import casadi.*

states   = SX.sym('states',nx,1);
controls = SX.sym('controls',nu,1);
params   = SX.sym('paras',np,1);

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
b=-m*l*cos(theta)*sin(theta)*omega^2+u*cos(theta)+(M+m)*g*sin(theta);
c=M+m-m*(cos(theta))^2;

x_dot=[v;omega;a/c;b/(l*c)];
xdot = SX.sym('xdot',nx,1);
impl_f = xdot - x_dot;
     
%% Objectives and constraints

% objectives
h = [p;theta;v;omega;u];

hN = h(1:nyN);

% general inequality constraints
general_con = [];
general_con_N = [];

%% NMPC sampling time [s]

Ts = 0.025; % simulation sample time
Ts_st = 0.025; % shooting interval time

%% build casadi function (don't touch)
h_fun=Function('h_fun', {states,controls,params}, {h},{'states','controls','params'},{'h'});
hN_fun=Function('hN_fun', {states,params}, {hN},{'states','params'},{'hN'});

path_con_fun=Function('path_con_fun', {states,controls,params}, {general_con},{'states','controls','params'},{'general_con'});
path_con_N_fun=Function('path_con_N_fun', {states,params}, {general_con_N},{'states','params'},{'general_con_N'});

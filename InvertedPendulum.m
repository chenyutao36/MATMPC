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
nc=2;
ncN=1;

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
% u=states(5);
% du=controls(1);

a=-m*l*sin(theta)*omega^2+m*g*cos(theta)*sin(theta)+u;
b=-m*l*cos(theta)*omega^2+u*cos(theta)+(M+m)*g*sin(theta);
c=M+m-m*(cos(theta))^2;

x_dot=[v;omega;a/c;b/(l*c)];
% x_dot=[v;omega;a/c;b/(l*c);du];

xdot = SX.sym('xdot',nx,1);
impl_f = xdot - x_dot;
     
%% Objectives and constraints

h = [p;theta;v;omega;u];
% h = [theta;du];

hN = h(1:nyN);

h_fun=Function('h_fun', {states,controls,params}, {h},{'states','controls','params'},{'h'});
hN_fun=Function('hN_fun', {states,params}, {hN},{'states','params'},{'hN'});

ineq=[p;u];
ineqN=p;
% ineq=u;
% ineqN=u;
ineq_fun=Function('ineq_fun', {states,controls,params}, {ineq},{'states','controls','params'},{'ineq'});
ineqN_fun=Function('ineqN_fun', {states,params}, {ineqN},{'states','params'},{'ineqN'});

lb_ineq=SX.sym('lb_ineq',length(ineq),1);
ub_ineq=SX.sym('ub_ineq',length(ineq),1);
lbN_ineq=SX.sym('lbN_ineq',length(ineqN),1);
ubN_ineq=SX.sym('ubN_ineq',length(ineqN),1);
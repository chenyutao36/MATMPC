%------------------------------------------%

% Model

%------------------------------------------%

%% Dimensions

nx=12;  % No. of states
nu=6;  % No. of controls
ny=6; % No. of outputs
nyN=6; % No. of outputs at the terminal point
np=0; % No. of model parameters
nc=0; % No. of general inequality constraints
ncN=0; % No. of general inequality constraints at the terminal point
nbx = 0; % No. of bounds on states
nbu = 6; % No. of bounds on controls

% state and control bounds
nbx_idx = 0;  % indexs of states which are bounded
nbu_idx = 1:6;  % indexs of controls which are bounded

%% create variables

import casadi.*

states   = SX.sym('states',nx,1);
controls = SX.sym('controls',nu,1);
params   = SX.sym('paras',np,1);

%% Dynamics

%Gravity
g=9.81;
%Distance between CoM and rotor i
l=0.315;
%Mass of Hexacopter
m=2.4;
%Constant
c_tau=1.7*10^-2;
%Inertia
Ix=0.04;
Iy=0.0422;
Iz=0.0827;
%Tilted angle
alpha=30*pi/180;

k1 = l*cos(alpha) + c_tau*sin(alpha);
k2 = -l*sin(alpha) + c_tau*cos(alpha);

phi=states(1);
theta=states(2);
psi=states(3);
p=states(4);  
q=states(5);
r=states(6);
u=states(7);
v=states(8);
w=states(9);
x=states(10);
y=states(11);
z=states(12);

f1=controls(1);
f2=controls(2);
f3=controls(3);
f4=controls(4);
f5=controls(5);
f6=controls(6);

x_dot = [ p+r*cos(phi)*tan(theta)+q*sin(phi)*tan(theta);...
          q*cos(phi) - r*sin(phi);...
          r*cos(phi)/cos(theta) + q*sin(phi)/cos(theta);...
          q*r*(Iy-Iz)/Ix + sqrt(3)/2*k1*(f2+f3-f5-f6)/Ix;...
          r*p*(Iz-Ix)/Iy + k1*(-f1-0.5*f2+0.5*f3+f4+0.5*f5-0.5*f6)/Iy;...
          q*p*(Ix-Iy)/Iz + k2*(f1-f2+f3-f4+f5-f6)/Iz;...
          g*sin(theta) + sqrt(3)*sin(alpha)/2/m*(-f2+f3-f5+f6) + r*v-q*w;...
          -g*sin(phi)*cos(theta) + sin(alpha)/m*(-f1+0.5*f2+0.5*f3-f4+0.5*f5+0.5*f6) + p*w-r*u;...
          -g*cos(theta)*cos(phi) + cos(alpha)/m*(f1+f2+f3+f4+f5+f6) + q*u-p*v;...
          w*(sin(phi)*sin(psi)+cos(phi)*cos(psi)*sin(theta)) - v*(cos(phi)*sin(psi)-cos(psi)*sin(phi)*sin(theta)) + u*cos(psi)*cos(theta);...
          v*(cos(phi)*cos(psi)+sin(phi)*sin(psi)*sin(theta)) - w*(sin(phi)*cos(psi)-sin(psi)*cos(phi)*sin(theta)) + u*sin(psi)*cos(theta);...
          w*cos(phi)*cos(theta)-u*sin(theta)+v*cos(theta)*sin(phi);
    ];


xdot = SX.sym('xdot',nx,1);
impl_f = xdot - x_dot;
     
%% Objectives and constraints

h = [x;y;z;phi;theta;psi];

hN = h(1:nyN);

% general inequality path constraints
general_con=[];  
general_con_N=[]; 

%% NMPC sampling time [s]

Ts = 0.01; % simulation sample time
Ts_st = 0.1; % shooting interval time

%% build casadi function (don't touch)

h_fun=Function('h_fun', {states,controls,params}, {h},{'states','controls','params'},{'h'});
hN_fun=Function('hN_fun', {states,params}, {hN},{'states','params'},{'hN'});

path_con_fun=Function('path_con_fun', {states,controls,params}, {general_con},{'states','controls','params'},{'general_con'});
path_con_N_fun=Function('path_con_N_fun', {states,params}, {general_con_N},{'states','params'},{'general_con_N'});



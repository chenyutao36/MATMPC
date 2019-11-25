%------------------------------------------%
% Tethered quadcopter for a safe and constrained maneuver

% from "Online Nonlinear Model Predictive Control for Tethered UAVs to
% Perform a Safe and Constrained Maneuver  ", E.Rossi, 2019

% typical configuration: 1) N=30, Ts=Ts_st=0.01, no shifting 

%------------------------------------------%


%% Dimensions

nx=6;  % No. of states
nu=4;  % No. of controls
nz=0;
ny=12; % No. of outputs
nyN=4; % No. of outputs at the terminal point
np=2; % No. of model parameters
nc=3;%0; % No. of general inequality constraints
ncN=1;%1; % No. of general inequality constraints
nbx = 2; % No. of bounds on states
nbu = 2; % No. of bounds on controls

% state and control bounds
nbx_idx = 5:6;  % indexs of states which are bounded
nbu_idx = 3:4;  % indexs of controls which are bounded

%% create variables

import casadi.*

states   = SX.sym('states',nx,1);
controls = SX.sym('controls',nu,1);
alg      = SX.sym('alg',nz,1);
params   = SX.sym('paras',np,1);
refs     = SX.sym('refs',ny,1);     % references of the first N stages
refN     = SX.sym('refs',nyN,1);    % reference of the last stage
Q        = SX.sym('Q',ny,1);        % weighting matrix of the first N stages
QN       = SX.sym('QN',nyN,1);      % weighting matrix of the last stage
aux      = SX.sym('aux',ny,1);      % auxilary variable
auxN     = SX.sym('auxN',nyN,1);    % auxilary variable


%% Dynamics

%Constant

g = 9.81; % gravity [m/s^2]

l = 1; % link length = 1[m]
mR = 1; % drone mass [Kg]
JR = 0.15; % drone inertia
a1 = -g/l;
a2 = 1/(mR*l);
a3 = 1/JR;
d = 0.2;

phi=states(1);
phi_dot=states(2);
theta=states(3);
theta_dot=states(4);
f1=states(5);%b*omega1^2
f2=states(6);%b*omega2^2
df1=controls(1);%b*omega1^2
df2=controls(2);%b*omega2^2
s1 = controls(3);
s2 = controls(4);
phi_ref = params(1);
theta_ref = params(2);

x_dot = [phi_dot; a1*cos(phi); theta_dot; 0] + ...
        [0                  0;
         a2*cos(phi+theta)  0;
         0                  0;
         0                  a3]  *  [f1; f2];%*   [1 1; -d d]

x_dot = [x_dot;
         df1;          
         df2];
     
z_fun = [];

xdot = SX.sym('xdot',nx,1);
impl_f = xdot - x_dot;
     
%% Objectives and constraints
phi_lim_vel = phi_ref + 20*pi/180;
phi_lim = phi_ref + 16*pi/180;
gamma2_vel = 180/pi/0.8;
gamma1 = 1;
gamma2 = 180/pi/3;

% inner objectives
h = [phi-phi_ref;phi_dot;theta-theta_ref;theta_dot;df1;df2;
     gamma1*(theta-theta_ref)*(1/(1+exp(-gamma2*(phi_lim-phi)))); % to arrive close to the ground with attitude = surface slope
     phi_dot*(1/(1+exp(-gamma2_vel*(phi_lim_vel-phi)))); % elev velocity goes to zero as the drone is closer to the surface
     theta_dot*(1/(1+exp(-gamma2_vel*(phi_lim_vel-phi)))); % attitude velocity goes to zero as the drone is closer to the surface
     1/(phi+theta-pi/2);...
     s1;s2]; 
hN = [phi-phi_ref;phi_dot;theta-theta_ref;theta_dot];

% outer objectives
obji = 0.5*(h-refs)'*diag(Q)*(h-refs);
objN = 0.5*(hN-refN)'*diag(QN)*(hN-refN);

obji_GGN = 0.5*(aux-refs)'*(aux-refs);
objN_GGN = 0.5*(auxN-refN)'*(auxN-refN);

% general inequality path constraints
general_con = [1/a2*phi_dot^2 + a1/a2*sin(phi) + sin(phi+theta)*(f1+f2); 
                d*sin(theta-pi/6)-l*sin(pi/6+phi)-s1;
                -d*sin(theta-pi/6)-l*sin(pi/6+phi)-s2];
general_con_N = [1/a2*phi_dot^2 + a1/a2*sin(phi) + sin(phi+theta)*(f1+f2)]; 

%% NMPC discretizing time length [s]

Ts_st = 0.01; % shooting interval time

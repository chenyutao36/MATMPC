%------------------------------------------%

% Model

%------------------------------------------%


%% Dimensions

nx=6;  % No. of states
nu=4;  % No. of controls
ny=12; % No. of outputs
nyN=4; % No. of outputs at the terminal point
np=2; % No. of model parameters
nc=3;%0; % No. of general inequality constraints
ncN=1;%1; % No. of general inequality constraints
nbx = 2; % No. of bounds on states
nbu = 2; % No. of bounds on controls

import casadi.*

states   = SX.sym('states',nx,1);
controls = SX.sym('controls',nu,1);
params   = SX.sym('paras',np,1);
refs     = SX.sym('refs',ny,1);
refN     = SX.sym('refs',nyN,1);
Q        = SX.sym('Q',ny,ny);
QN       = SX.sym('QN',nyN,nyN);


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

xdot = SX.sym('xdot',nx,1);
impl_f = xdot - x_dot;
     
%% Objectives and constraints
% alpha = pi/6;
% phi_ref = -alpha;
% theta_ref = alpha;
phi_lim_vel = phi_ref + 20*pi/180;
phi_lim = phi_ref + 16*pi/180;
gamma2_vel = 180/pi/0.8;
gamma1 = 1;
gamma2 = 180/pi/3;


h = [phi-phi_ref;phi_dot;theta-theta_ref;theta_dot;df1;df2;
     gamma1*(theta-theta_ref)*(1/(1+exp(-gamma2*(phi_lim-phi)))); % to arrive close to the ground with attitude = surface slope
     phi_dot*(1/(1+exp(-gamma2_vel*(phi_lim_vel-phi)))); % elev velocity goes to zero as the drone is closer to the surface
     theta_dot*(1/(1+exp(-gamma2_vel*(phi_lim_vel-phi)))); % attitude velocity goes to zero as the drone is closer to the surface
     1/(phi+theta-pi/2);...
     s1;s2]; 

hN = [phi-phi_ref;phi_dot;theta-theta_ref;theta_dot];

h_fun=Function('h_fun', {states,controls,params}, {h},{'states','controls','params'},{'h'});
hN_fun=Function('hN_fun', {states,params}, {hN},{'states','params'},{'hN'});

% general inequality path constraints
general_con = [1/a2*phi_dot^2 + a1/a2*sin(phi) + sin(phi+theta)*(f1+f2); 
                d*sin(theta-pi/6)-l*sin(pi/6+phi)-s1;
                -d*sin(theta-pi/6)-l*sin(pi/6+phi)-s2];
general_con_N = [1/a2*phi_dot^2 + a1/a2*sin(phi) + sin(phi+theta)*(f1+f2)]; 
%     d*sin(theta-pi/6)-l*sin(pi/6+phi)+s1;
%     -d*sin(theta-pi/6)-l*sin(pi/6+phi)+s2];

% state and control bounds
nbx_idx = 5:6;  % indexs of states which are bounded
nbu_idx = 3:4;  % indexs of controls which are bounded
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

Ts = 0.01; % simulation sample time
Ts_st = 0.01; % shooting interval time
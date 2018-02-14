%------------------------------------------%

% Model TiltHex

%------------------------------------------%


%% Dimensions

nx=18;  % No. of states
nu=6;  % No. of controls
ny=12; % No. of outputs
nyN=6; % No. of outputs at the terminal point
np=0; % No. of model parameters
nc=0; % No. of general inequality constraints
ncN=0; % No. of general inequality constraints at the terminal point
nbx = 6; % No. of bounds on states
nbu = 6; % No. of bounds on controls

import casadi.*

states   = SX.sym('states',nx,1);
controls = SX.sym('controls',nu,1);
params   = SX.sym('paras',np,1);
refs     = SX.sym('refs',ny,1);
refN     = SX.sym('refs',nyN,1);
Q        = SX.sym('Q',ny,ny);
QN       = SX.sym('QN',nyN,nyN);


%% Dynamics

%Gravity
g=9.81;
%Distance between CoM and rotor i
l=0.4;
%l=0.315;
%Mass of TiltHex
m=1.8;
% m=2.4;
%Constant
c_tau=1.7*10^-2;
%Inertia
% Ix=0.04;
% Iy=0.0422;
% Iz=0.0827;
%Inertia
% Ix=11.5e-6;
% Iy=11.4e-6;
% Iz=19.4e-6;
Ix=11549*10^(-6);
Iy=11368*10^(-6);
Iz=19444*10^(-6);
%Tilted angle
alpha=35*pi/180;
beta=-10*pi/180;

%Costants for the thrust
fx1 = 0.5*sin(beta)-sqrt(3)*0.5*cos(beta)*sin(alpha);
fx2 = 0.5*sin(beta)+sqrt(3)*0.5*cos(beta)*sin(alpha);
fy1 = sqrt(3)*0.5*sin(beta)+0.5*cos(beta)*sin(alpha);
fy2 = -sqrt(3)*0.5*sin(beta)+0.5*cos(beta)*sin(alpha);
taux1 = -0.5*c_tau*sin(beta)+sqrt(3)*0.5*(l*cos(alpha)+c_tau*sin(alpha))*cos(beta);
taux2 = -0.5*c_tau*sin(beta)-sqrt(3)*0.5*(l*cos(alpha)+c_tau*sin(alpha))*cos(beta);
tauy1 = (l*cos(alpha)+c_tau*sin(alpha))*cos(beta);
tauy2 = sqrt(3)*0.5*c_tau*sin(beta)+0.5*(l*cos(alpha)+c_tau*sin(alpha))*cos(beta);
tauy3 = -sqrt(3)*0.5*c_tau*sin(beta)+0.5*(l*cos(alpha)+c_tau*sin(alpha))*cos(beta);
tauz1 = (l*sin(alpha)-c_tau*cos(alpha))*cos(beta);

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
f1=states(13);
f2=states(14);
f3=states(15);
f4=states(16);
f5=states(17);
f6=states(18);

df1=controls(1);
df2=controls(2);
df3=controls(3);
df4=controls(4);
df5=controls(5);
df6=controls(6);

x_dot = [ p+r*cos(phi)*tan(theta)+q*sin(phi)*tan(theta);...
          q*cos(phi) - r*sin(phi);...
          r*cos(phi)/cos(theta) + q*sin(phi)/cos(theta);...
          q*r*(Iy-Iz)/Ix + (c_tau*sin(beta)*(f1+f4)+taux1*(f2+f3)+taux2*(f5+f6))/Ix;...
          r*p*(Iz-Ix)/Iy + (tauy1*(-f1+f4)+tauy2*(-f2+f3)+tauy3*(f5-f6))/Iy;...
          q*p*(Ix-Iy)/Iz + tauz1*(-f1+f2-f3+f4-f5+f6)/Iz;...
          g*sin(theta) + (sin(beta)*(f1-f4)+fx1*(f2-f3)+fx2*(f6-f5))/m + r*v-q*w;...
          -g*sin(phi)*cos(theta) + (cos(beta)*sin(alpha)*(-f1-f4)+fy1*(f2+f3)+fy2*(f5+f6))/m + p*w-r*u;...
          -g*cos(theta)*cos(phi) + cos(beta)*cos(alpha)*(f1+f2+f3+f4+f5+f6)/m + q*u-p*v;...
          w*(sin(phi)*sin(psi)+cos(phi)*cos(psi)*sin(theta)) - v*(cos(phi)*sin(psi)-cos(psi)*sin(phi)*sin(theta)) + u*cos(psi)*cos(theta);...
          v*(cos(phi)*cos(psi)+sin(phi)*sin(psi)*sin(theta)) - w*(sin(phi)*cos(psi)-sin(psi)*cos(phi)*sin(theta)) + u*sin(psi)*cos(theta);...
          w*cos(phi)*cos(theta)-u*sin(theta)+v*cos(theta)*sin(phi);
          df1;...
          df2;...
          df3;...
          df4;...
          df5;...
          df6
    ];


xdot = SX.sym('xdot',nx,1);
impl_f = xdot - x_dot;
     
%% Objectives and constraints

h = [x;y;z;phi;theta;psi;controls];

hN = h(1:nyN);

h_fun=Function('h_fun', {states,controls,params}, {h},{'states','controls','params'},{'h'});
hN_fun=Function('hN_fun', {states,params}, {hN},{'states','params'},{'hN'});

% general inequality path constraints
general_con=[];  
general_con_N=[]; 

% state and control bounds
nbx_idx = 13:18;  % indexs of states which are bounded
nbu_idx = 1:6;  % indexs of controls which are bounded
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

Ts = 0.002; % simulation sample time
Ts_st = 0.1; % shooting interval time
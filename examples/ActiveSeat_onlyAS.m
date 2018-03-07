%% find your path to the original active seat model files

addpath(genpath('C:\Users\enrico\Documents\MATLAB\active seat(original)\nonlinear'));

cd('C:\Users\enrico\Documents\MATLAB\active seat(original)\nonlinear');

%% initialization parameters
m = 67;
sigma_0 = 10^4;
vs = 0.005;
Fs = 45;
Fc = 30;
g = 9.81;
alpha = 10;
MM = 50;
k1 = 12000;
k2 = 1000;
c1 = 200;
c2 = 2000;

%% Dimensions

nx=8;   % No. of states
nu=4;   % No. of controls
ny=7;   % No. of outputs
nyN=1;  % No. of outputs at the terminal point
np=0;   % No. of model parameters
nc=0;   % No. of general constraints
ncN=0;  % No. of general constraints at the terminal point
nbx = 0;% nbx = 6; % No. of bounds on states
nbu = 0;% No. of bounds on controls


import casadi.*

states   = SX.sym('states',nx,1);
controls = SX.sym('controls',nu,1);
params   = SX.sym('paras',np,1);
refs     = SX.sym('refs',ny,1);
refN     = SX.sym('refs',nyN,1);
Q        = SX.sym('Q',ny,ny);
QN       = SX.sym('QN',nyN,nyN);


%% Active Seat Dynamics

accX=states(1);
roll=states(2);
accY=states(3);
prY1=states(4); 
prY2=states(5); 
prY3=states(6); 
y_press=states(7);
pressY=states(8);

daccX=controls(1);
droll=controls(2);
daccY=controls(3);
dpressY=controls(4);

tmp1= (sqrt(prY2^2)*prY3) ;
tmp2= m*accX*cos(pi/180*alpha)+MM*g*sin(pi/180*alpha) ;
tmp3= 1/(pi)*atan(tmp2)+0.6 ;

x_dot=[ daccX;...
        droll;...
        daccY;...
        prY2;...
       -(c1*(10*prY1)^2+c2)/m*prY2-(k1*(10*prY1)^2+k2)*prY1/m+accY+g*roll-sigma_0*prY3/m; ...
       prY2-tmp1/((Fc*tmp3+((Fs-Fc)*tmp3*exp(-(prY2/vs)^2)))/sigma_0);...               
       200*k1*prY1^2*prY2+(100*k1*prY1^2+k2)*prY2+dpressY;...
       dpressY];
 
xdot = SX.sym('xdot',nx,1);
impl_f = xdot - x_dot;
     
%% Objectives and constraints

% objectives
h = [y_press; controls(1:3);pressY; accX; accY ]; 
hN=pressY;

h_fun=Function('h_fun', {states,controls,params}, {h},{'states','controls','params'},{'h'});
hN_fun=Function('hN_fun', {states,params}, {hN},{'states','params'},{'hN'});

% general inequality constraints
general_con = [];
general_con_N = [];

% state and control bounds
nbx_idx = 0;  % nbx_idx = [2]; %nbx_idx = [5,6,14,15,21,26]; % indexs of states which are bounded
nbu_idx = 0;    % indexs of controls which are bounded
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

Ts = 0.005; % simulation sample time
Ts_st = 0.005; % shooting interval time

%%

cd('C:\Users\enrico\Documents\MATLAB\GITLAB\NMPCtool');

clc;
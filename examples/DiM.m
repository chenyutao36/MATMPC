%------------------------------------------%
% Motion cueing model for a  six-dof dynamic driving simulator

% modified from "A Nonlinear, MPC-Based Motion Cueing Algorithm for a 
% High-Performance, Nine-DOF Dynamic Simulator Platform", Bruschetta, 2016 

%------------------------------------------%

%% Filters
cutoff_frequency=2;  % cuting-off frequence=2 Hz;
sampling_frequency=100;

[A_hp,B_hp,C_hp,D_hp]=butter(1,cutoff_frequency/sampling_frequency*2,'high');
[A_lp,B_lp,C_lp,D_lp]=butter(1,cutoff_frequency/sampling_frequency*2);

hp_c=d2c(ss(A_hp,B_hp,C_hp,D_hp,1/sampling_frequency),'zoh');
lp_c=d2c(ss(A_lp,B_lp,C_lp,D_lp,1/sampling_frequency),'zoh');

[A_hp_c,B_hp_c,C_hp_c,D_hp_c] = ssdata(hp_c);
[A_lp_c,B_lp_c,C_lp_c,D_lp_c] = ssdata(lp_c);

n=size(A_hp_c,1);

%% Dimensions

nx=30+6*n;  % No. of states
nu=6;  % No. of controls
ny=30; % No. of outputs
nyN=24; % No. of outputs at the terminal point
np=0; % No. of model parameters
nc=6; % No. of general inequality constraints
ncN=6; % No. of general inequality constraints at the terminal point
nbx = 0; % No. of bounds on states
nbu = 0; % No. of bounds on controls

% state and control bounds
nbx_idx = 0;  % indexs of states which are bounded
nbu_idx = 0;  % indexs of controls which are bounded

%% create variables

import casadi.*

states   = SX.sym('states',nx,1);   % state variables
controls = SX.sym('controls',nu,1); % control inputs
params   = SX.sym('paras',np,1);    % parameters

%% Vestibular model

[A_S, B_S, C_S, D_S]=tf2ss([5.73*80, 0 0],[5.73*80, 5.73+80, 1]);  % SCC state-space model for angular rotation
[A_O, B_O, C_O, D_O]=tf2ss([0.4*10, 0.4*1],[0.016*5, 0.016+5, 1]); % OTH state-space model for linear motion

A_hat=blkdiag(A_O,A_O,A_O,A_S,A_S,A_S);  % a second order model for each of 3 DOFs
B_hat=blkdiag(B_O,B_O,B_O,B_S,B_S,B_S);
C_hat=blkdiag(C_O,C_O,C_O,C_S,C_S,C_S);
D_hat=blkdiag(D_O,D_O,D_O,D_S,D_S,D_S);

%% Dynamics

g=9.81;

x_OTH_1=states(1);  
x_OTH_2=states(2);  
y_OTH_1=states(3);  
y_OTH_2=states(4);  
z_OTH_1=states(5);  
z_OTH_2=states(6);  

psi_SCC_1=states(7); 
psi_SCC_2=states(8); 
theta_SCC_1=states(9); 
theta_SCC_2=states(10); 
phi_SCC_1=states(11); 
phi_SCC_2=states(12); 

x_hex_tri=states(13);
y_hex_tri=states(14);
z=states(15);

vx_hex_tri=states(16);
vy_hex_tri=states(17);
vz=states(18);

x_tri=states(19);
y_tri=states(20);

vx_tri=states(21);
vy_tri=states(22);

phi_hex=states(23);
phi_hex_d=states(24);

theta_hex=states(25);
theta_hex_d=states(26);

psi_hex=states(27);
psi_hex_d=states(28);

phi_tri=states(29);
phi_tri_d=states(30);

ax=controls(1);
ay=controls(2);
az=controls(3);
phi_dd=controls(4);
theta_dd=controls(5);
psi_dd=controls(6);

% rotation matrix
R_tri=[cos(phi_tri) -sin(phi_tri) 0;...
       sin(phi_tri) cos(phi_tri)  0;...
       0            0             1];

R_hex=[cos(phi_hex) -sin(phi_hex) 0;...
       sin(phi_hex) cos(phi_hex)  0;...
       0            0             1]*...
      [cos(theta_hex)  0  sin(theta_hex);...
       0               1  0;...
       -sin(theta_hex) 0  cos(theta_hex)]*...
      [1  0             0;...
       0  cos(psi_hex)  -sin(psi_hex);...
       0  sin(psi_hex)  cos(psi_hex)];
   
R_tri_hex=R_tri*R_hex;

R_tri_hex_inv=R_tri_hex';

ax_hex_tri=C_hp_c*states(30+0*n+1:30+1*n)+D_hp_c*ax;
ay_hex_tri=C_hp_c*states(30+1*n+1:30+2*n)+D_hp_c*ay;

ax_tri=C_lp_c*states(30+2*n+1:30+3*n)+D_lp_c*ax;
ay_tri=C_lp_c*states(30+3*n+1:30+4*n)+D_lp_c*ay;

phi_hex_dd=C_hp_c*states(30+4*n+1:30+5*n)+D_hp_c*phi_dd;
phi_tri_dd=C_lp_c*states(30+5*n+1:30+6*n)+D_lp_c*phi_dd;

% W=[1 0 -sin(theta_hex);...
%    0 cos(psi_hex) sin(psi_hex)*cos(theta_hex);
%    0 -sin(psi_hex) cos(psi_hex)*cos(theta_hex)];
W=eye(3);

omega_tri=[0;0;phi_tri_d];
omega_hex=W*[psi_hex_d;theta_hex_d;phi_hex_d];
omega_tri_hex=W*omega_tri+omega_hex;

alpha_tri=[0;0;phi_tri_dd];
% alpha_tri_hex=[psi_dd-theta_hex_d*(phi_hex_dd+phi_tri_dd);...
%                theta_dd+psi_hex_d*(phi_hex_dd+phi_tri_dd);...
%                (phi_hex_dd+phi_tri_dd)-psi_hex_d*theta_hex_d];

a_d=[ax_tri;...
     ay_tri;...
     0]+...
     R_tri*[ax_hex_tri;...
            ay_hex_tri;...
            az]+...
     2*cross(omega_tri,R_tri*[vx_hex_tri;vy_hex_tri;vz]); ...
%      + cross( alpha_tri_hex, R_tri_hex*[0.5;0.3;0.8] );
 %      + cross( alpha_tri, R_tri*[x_hex_tri;y_hex_tri;z] )...
%      + cross( omega_tri_hex, cross(omega_tri_hex,R_tri_hex*[0.5;0.3;1]) )...
%      + cross( omega_tri, cross(omega_tri,R_tri*[x_hex_tri;y_hex_tri;z]) )...

u_VEST=[R_tri_hex_inv*a_d+[R_hex(1:2,:);0 0 0]*[0;0;-g];...
        omega_tri_hex];

 
y_VEST=C_hat*states(1:12)+D_hat*u_VEST;

x_VEST=A_hat*states(1:12)+B_hat*u_VEST;

% dynamics
x_dot=[x_VEST;vx_hex_tri;...
             vy_hex_tri;...
             vz;...
             ax_hex_tri;...
             ay_hex_tri;... 
             az;...
             vx_tri;...
             vy_tri;...
             ax_tri;...
             ay_tri;...
             phi_hex_d;...
             phi_hex_dd;...
             theta_hex_d;...
             theta_dd;...
             psi_hex_d;...
             psi_dd;...
             phi_tri_d;...
             phi_tri_dd;...
             A_hp_c*states(31:30+n)+B_hp_c*ax;...
             A_hp_c*states(30+n+1:30+2*n)+B_hp_c*ay;...
             A_lp_c*states(30+2*n+1:30+3*n)+B_lp_c*ax;...
             A_lp_c*states(30+3*n+1:30+4*n)+B_lp_c*ay;...
             A_hp_c*states(30+4*n+1:30+5*n)+B_hp_c*phi_dd;...
             A_lp_c*states(30+5*n+1:30+6*n)+B_lp_c*phi_dd;...
         ];
 
xdot = SX.sym('xdot',nx,1);
impl_f = xdot - x_dot;

%% Objectives and constraints

h = [y_VEST;...
     x_hex_tri; y_hex_tri; z;...
     vx_hex_tri; vy_hex_tri; vz;...
     x_tri; y_tri;...
     vx_tri; vy_tri;...
     phi_hex; theta_hex; psi_hex;...
     phi_tri;...
     phi_hex_d; theta_hex_d;psi_hex_d;...
     phi_tri_d;...
     ax;
     ay;
     az;...
     phi_dd;...
     theta_dd;...
     psi_dd...
    ];

hN = h(1:ny-nu);

% general inequality path constraints
q1 = ((x_hex_tri + (2126791213464199*cos(-theta_hex)*cos(-phi_hex))/9007199254740992 - (5935510252496305*cos(-psi_hex)*sin(-phi_hex))/9007199254740992 - (5935510252496305*cos(-phi_hex)*sin(-theta_hex)*sin(-psi_hex))/9007199254740992 + 8346791275420823/18014398509481984)^2 + (y_hex_tri + (5935510252496305*cos(-psi_hex)*cos(-phi_hex))/9007199254740992 + (2126791213464199*cos(-theta_hex)*sin(-phi_hex))/9007199254740992 - (5935510252496305*sin(-theta_hex)*sin(-psi_hex)*sin(-phi_hex))/9007199254740992 - 592520818642875/562949953421312)^2 + (z + (2126791213464199*sin(-theta_hex))/9007199254740992 + (5935510252496305*cos(-theta_hex)*sin(-psi_hex))/9007199254740992 + 1811/2000)^2)^(1/2);
q2 = (((cos(-psi_hex)*sin(-phi_hex))/8 - x_hex_tri + (cos(-phi_hex)*sin(-theta_hex)*sin(-psi_hex))/8 + (759^(1/2)*cos(-theta_hex)*cos(-phi_hex))/40 - 765438935074863/1125899906842624)^2 + (z + (cos(-theta_hex)*sin(-psi_hex))/8 - (759^(1/2)*sin(-theta_hex))/40 + 1811/2000)^2 + ((sin(-theta_hex)*sin(-psi_hex)*sin(-phi_hex))/8 - (cos(-psi_hex)*cos(-phi_hex))/8 - y_hex_tri + (759^(1/2)*cos(-theta_hex)*sin(-phi_hex))/40 + 522152074465211/562949953421312)^2)^(1/2);
q3 = ((x_hex_tri + (cos(-psi_hex)*sin(-phi_hex))/8 + (cos(-phi_hex)*sin(-theta_hex)*sin(-psi_hex))/8 - (759^(1/2)*cos(-theta_hex)*cos(-phi_hex))/40 + 6123511480598909/9007199254740992)^2 + (y_hex_tri - (cos(-psi_hex)*cos(-phi_hex))/8 + (sin(-theta_hex)*sin(-psi_hex)*sin(-phi_hex))/8 - (759^(1/2)*cos(-theta_hex)*sin(-phi_hex))/40 + 8354433191443371/9007199254740992)^2 + (z - (cos(-theta_hex)*sin(-psi_hex))/8 - (759^(1/2)*sin(-theta_hex))/40 + 1811/2000)^2)^(1/2);
q4 = ((z + (8507164853856779*sin(-theta_hex))/36028797018963968 - (5935510252496307*cos(-theta_hex)*sin(-psi_hex))/9007199254740992 + 1811/2000)^2 + (x_hex_tri + (8507164853856779*cos(-theta_hex)*cos(-phi_hex))/36028797018963968 + (5935510252496307*cos(-psi_hex)*sin(-phi_hex))/9007199254740992 + (5935510252496307*cos(-phi_hex)*sin(-theta_hex)*sin(-psi_hex))/9007199254740992 + 260837227356901/562949953421312)^2 + (y_hex_tri - (5935510252496307*cos(-psi_hex)*cos(-phi_hex))/9007199254740992 + (8507164853856779*cos(-theta_hex)*sin(-phi_hex))/36028797018963968 + (5935510252496307*sin(-theta_hex)*sin(-psi_hex)*sin(-phi_hex))/9007199254740992 + 4740166549142999/4503599627370496)^2)^(1/2);
q5 = ((x_hex_tri + (1019226764088171*cos(-theta_hex)*cos(-phi_hex))/2251799813685248 + (4809610345653685*cos(-psi_hex)*sin(-phi_hex))/9007199254740992 - 2091^(1/2)/40 + (4809610345653685*cos(-phi_hex)*sin(-theta_hex)*sin(-psi_hex))/9007199254740992)^2 + (y_hex_tri - (4809610345653685*cos(-psi_hex)*cos(-phi_hex))/9007199254740992 + (1019226764088171*cos(-theta_hex)*sin(-phi_hex))/2251799813685248 + (4809610345653685*sin(-theta_hex)*sin(-psi_hex)*sin(-phi_hex))/9007199254740992 + 1/8)^2 + (z + (1019226764088171*sin(-theta_hex))/2251799813685248 - (4809610345653685*cos(-theta_hex)*sin(-psi_hex))/9007199254740992 + 1811/2000)^2)^(1/2);
q6 = ((y_hex_tri + (300600646603355*cos(-psi_hex)*cos(-phi_hex))/562949953421312 + (8153814112705379*cos(-theta_hex)*sin(-phi_hex))/18014398509481984 - (300600646603355*sin(-theta_hex)*sin(-psi_hex)*sin(-phi_hex))/562949953421312 - 1/8)^2 + (z + (8153814112705379*sin(-theta_hex))/18014398509481984 + (300600646603355*cos(-theta_hex)*sin(-psi_hex))/562949953421312 + 1811/2000)^2 + ((300600646603355*cos(-psi_hex)*sin(-phi_hex))/562949953421312 - (8153814112705379*cos(-theta_hex)*cos(-phi_hex))/18014398509481984 - x_hex_tri + 2091^(1/2)/40 + (300600646603355*cos(-phi_hex)*sin(-theta_hex)*sin(-psi_hex))/562949953421312)^2)^(1/2);

general_con=[q1;q2;q3;q4;q5;q6];  
general_con_N=[q1;q2;q3;q4;q5;q6]; 

%% NMPC sampling time [s]

Ts = 0.01; % simulation sample time
Ts_st = 0.01; % shooting interval time

%% build casadi function (don't touch)

h_fun=Function('h_fun', {states,controls,params}, {h},{'states','controls','params'},{'h'});
hN_fun=Function('hN_fun', {states,params}, {hN},{'states','params'},{'hN'});

path_con_fun=Function('path_con_fun', {states,controls,params}, {general_con},{'states','controls','params'},{'general_con'});
path_con_N_fun=Function('path_con_N_fun', {states,params}, {general_con_N},{'states','params'},{'general_con_N'});


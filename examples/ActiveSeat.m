%% find your path to the original active seat model files

addpath(genpath('e:/study/NNID/Active_seat_belt/Nonlinear'));

cd('e:/study/NNID/Active_seat_belt/Nonlinear');

%% read data

filename = mfilename;
filename_marg = [filename '.xlsx'];
filen = 'Calabogie';
sdata = resreader('C_Segment_Calabogie_Drive_Close_200Hz.res');

GenerateInitialValueForCoder_test;

global Var_xy_ Var_yx_ Var_za_ Var_zv_ status_xy status_yx status_za status_zv;
Var_xy_ = Var_xy;
Var_yx_ = Var_yx;
Var_za_ = Var_za;
Var_zv_ = Var_zv;
status_xy = 0; 
status_yx = 0;
status_za = 0;
status_zv = 0;

time = sdata.time_TIME;
chassis_accelerations_longitudinal = sdata.chassis_accelerations_longitudinal; 
chassis_accelerations_lateral = sdata.chassis_accelerations_lateral;
chassis_accelerations_vertical = sdata.chassis_accelerations_vertical;
chassis_velocities_yaw = sdata.chassis_velocities_yaw;
chassis_displacements_pitch = sdata.chassis_displacements_pitch;
chassis_displacements_roll = sdata.chassis_displacements_roll;
chassis_displacements_roll_wrt_road = sdata.chassis_displacements_roll_wrt_road;
chassis_displacements_pitch_wrt_road = sdata.chassis_displacements_pitch_wrt_road;
chassis_displacements_vertical_wrt_road = sdata.chassis_displacements_vertical_wrt_road;

IN1_XY = interp1(time,chassis_accelerations_longitudinal,time(1):Var_xy.Ts:time(end),'spline');
IN1_YX = -interp1(time,chassis_accelerations_lateral,time(1):Var_xy.Ts:time(end),'spline');
IN1_ZV = interp1(time,chassis_accelerations_vertical,time(1):Var_xy.Ts:time(end),'spline');

IN1_ZA = interp1(time,chassis_velocities_yaw,time(1):Var_xy.Ts:time(end),'spline');
IN2_XY = interp1(time,chassis_displacements_pitch(1:end),time(1):Var_xy.Ts:time(end),'spline');
IN2_YX = -interp1(time,chassis_displacements_roll(1:end),time(1):Var_xy.Ts:time(end),'spline');

A_vestXY = Var_xy_.A_vest_XY;
B_vestXY = Var_xy_.B_vest_XY;
C_vestXY = Var_xy_.C_vest_XY;
D_vestXY = Var_xy_.D_vest_XY;

A_vestYX = Var_yx_.A_vest_YX;
B_vestYX = Var_yx_.B_vest_YX;
C_vestYX = Var_yx_.C_vest_YX;
D_vestYX = Var_yx_.D_vest_YX;

A_vestZvert = Var_zv_.A_vest_Zvert;
B_vestZvert = Var_zv_.B_vest_Zvert;
C_vestZvert = Var_zv_.C_vest_Zvert;
D_vestZvert = Var_zv_.D_vest_Zvert;

A_vestZang = Var_za_.A_vest_Zang;
B_vestZang = Var_za_.B_vest_Zang;
C_vestZang = Var_za_.C_vest_Zang;
D_vestZang = Var_za_.D_vest_Zang;

Avxy = [A_vestXY B_vestXY; zeros(2,7) zeros(2,2)];
Bvxy = [zeros(7,2); eye(2,2)];
Cvxy = [C_vestXY D_vestXY];
Dvxy = D_vestXY*0;

Avyx = [A_vestYX B_vestYX; zeros(2,7) zeros(2,2)];
Bvyx = [zeros(7,2); eye(2,2)];
Cvyx = [C_vestYX D_vestYX];
Dvyx = D_vestYX*0;

Avzv = [A_vestZvert B_vestZvert; zeros(1,4) zeros(1,1)];
Bvzv = [zeros(4,1); eye(1,1)];
Cvzv = [C_vestZvert D_vestZvert];
Dvzv = D_vestZvert*0;

Avza = [A_vestZang B_vestZang; zeros(1,3) zeros(1,1)];
Bvza = [zeros(3,1); eye(1,1)];
Cvza = [C_vestZang D_vestZang];
Dvza = D_vestZang*0;

A_p = blkdiag(Avxy,Avyx,Avzv,Avza); 
B_p = blkdiag(Bvxy,Bvyx,Bvzv,Bvza);
C_p = blkdiag(Cvxy,Cvyx,Cvzv,Cvza); 
D_p = blkdiag(Dvxy,Dvyx,Dvzv,Dvza);

A = [A_p, zeros(length(A_p),3); zeros(3,length(A_p)+3)];
B = [B_p, zeros(length(B_p),1); zeros(3,7)];
C = [C_p, zeros(size(C_p,1),3)]; %; zeros(3,27) zeros(3,6)];
D = [D_p, zeros(size(D_p,1),1)];

A=[A,zeros(30,2)];

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

nx=32;  % No. of states
nu=7;  % No. of controls
ny=32; % No. of outputs
nyN=18; % No. of outputs at the terminal point
np=0; % No. of model parameters
nc=0; % No. of general constraints
ncN=0; % No. of general constraints at the terminal point
nbx = 6; % No. of bounds on states
nbu = 0; % No. of bounds on controls


import casadi.*

states   = SX.sym('states',nx,1);
controls = SX.sym('controls',nu,1);
params   = SX.sym('paras',np,1);
refs     = SX.sym('refs',ny,1);
refN     = SX.sym('refs',nyN,1);
Q        = SX.sym('Q',ny,ny);
QN       = SX.sym('QN',nyN,nyN);


%% Dynamics

x1SccX=states(1);
x2SccX=states(2);  
x1OtoX=states(3);  
x2OtoX=states(4);  
pitch=states(5);  
pX=states(6);  
vX=states(7); 
accX=states(8); 
velPitch=states(9); 
x1SccY=states(10); 
x2SccY=states(11); 
x1OtoY=states(12); 
x2OtoY=states(13);
roll=states(14);
pY=states(15);
vY=states(16);
accY=states(17);
velRoll=states(18);
x1OtoZ=states(19);
x2OtoZ=states(20);
pZ=states(21);
vZ=states(22);
accZ=states(23);
x1SccZ=states(24);
x2SccZ=states(25);
yaw=states(26);
velYaw=states(27);
prY1=states(28);
prY2=states(29);
prY3=states(30);
y_press=states(31);
pressY=states(32);

daccX=controls(1);
dvelPitch=controls(2);
daccY=controls(3);
dvelRoll=controls(4);
daccZ=controls(5);
dvelYaw=controls(6);
dpressY=controls(7);

tmp1= (sqrt(prY2^2)*prY3) ;
tmp2= m*accX*cos(pi/180*alpha)+MM*g*sin(pi/180*alpha) ;
tmp3= 1/(pi)*atan(tmp2)+0.6 ;

x_dot=[A(1:27,:)*states+B(1:27,:)*controls;...
       prY2;...
       -(c1*(10*prY1)^2+c2)/m*prY2-(k1*(10*prY1)^2+k2)*prY1/m+accY+g*roll-sigma_0*prY3/m; ...
       prY2-tmp1/((Fc*tmp3+((Fs-Fc)*tmp3*exp(-(prY2/vs)^2)))/sigma_0);...               
       200*k1*prY1^2*prY2+(100*k1*prY1^2+k2)*prY2+dpressY;...
       dpressY];
 
xdot = SX.sym('xdot',nx,1);
impl_f = xdot - x_dot;
     
%% Objectives and constraints

C_bar = C;
C_bar = [C_bar,zeros(18,2)];

% objectives
h = [C_bar*states + D*controls; y_press; controls(1:6);pressY; accX; velPitch; accY; velRoll; accZ; velYaw];

hN = C_bar*states; 

h_fun=Function('h_fun', {states,controls,params}, {h},{'states','controls','params'},{'h'});
hN_fun=Function('hN_fun', {states,params}, {hN},{'states','params'},{'hN'});

% general inequality constraints
general_con = [];
general_con_N = [];

% state and control bounds
nbx_idx = [5,6,14,15,21,26]; % indexs of states which are bounded
nbu_idx = 0; % indexs of controls which are bounded
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

%% save your data in the path of your MATMPC

data_AS.y_lim_xy_pa = y_lim_xy_pa;
data_AS.y_lim_xy_pl = y_lim_xy_pl;
data_AS.y_lim_yx_pa = y_lim_yx_pa;
data_AS.y_lim_yx_pl = y_lim_yx_pl;
data_AS.y_lim_zv_pv = y_lim_zv_pv;
data_AS.y_lim_za_pa = y_lim_za_pa;

save('e:/study/NNID/MATMPC/data/ActiveSeat/data_AS.mat','data_AS');

cd('e:/study/NNID/MATMPC');

clc;
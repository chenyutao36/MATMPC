%------------------------------------------%
% Chain of masses connected by linear springs

% from "Efficient direct multiple shooting for nonlinear model predictive
% control on long horizons", Kirches, 2012

% typical configuration: 1) N=40(80,160), Ts=Ts_st=0.2, no shifting 

%------------------------------------------%

%% Dimensions

n=5; 
nx=n*3+(n-1)*3;
nu=3;
nz=0;
np=0;
ny=3*(n-1)+3+3;
nyN=3*(n-1)+3;
nc=0;
ncN=0;
nbx=0;
nbu=3;

% state and control bounds
nbx_idx = 0;  % indexs of states which are bounded
nbu_idx = 1:3;  % indexs of controls which are bounded

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

k=0.1;
lr=0.55;
m=0.45;
g=9.81;

p0x=0;p0y=0;p0z=0;
xend=1;yend=0;zend=0;

px=states(1:n,1);
py=states(n+1:2*n,1);
pz=states(2*n+1:3*n,1);
vx=states(3*n+1:4*n-1,1);
vy=states(4*n:5*n-2,1);
vz=states(5*n-1:6*n-3,1);
ux=controls(1);
uy=controls(2);
uz=controls(3);

dist=SX(n,1);
Fx=SX(n,1);
Fy=SX(n,1);
Fz=SX(n,1);
ax=SX(n-1,1);
ay=SX(n-1,1);
az=SX(n-1,1);

dist(1)= ((px(1)-p0x)^2+(py(1)-p0y)^2+(pz(1)-p0z)^2)^0.5 ;
Fx(1)=(px(1)-p0x)*k*(n-lr/dist(1));
Fy(1)=(py(1)-p0y)*k*(n-lr/dist(1));
Fz(1)=(pz(1)-p0z)*k*(n-lr/dist(1));
for i=2:n
    dist(i)= ((px(i)-px(i-1))^2+(py(i)-py(i-1))^2+(pz(i)-pz(i-1))^2)^0.5 ;
    Fx(i)= (px(i)-px(i-1))*k*(n-lr/dist(i)) ;
    Fy(i)= (py(i)-py(i-1))*k*(n-lr/dist(i)) ;
    Fz(i)= (pz(i)-pz(i-1))*k*(n-lr/dist(i)) ;
    ax(i-1,1)= (Fx(i)-Fx(i-1))*n/m ;
    ay(i-1,1)= (Fy(i)-Fy(i-1))*n/m ;
    az(i-1,1)= (Fz(i)-Fz(i-1))*n/m-g ;
end

x_dot=[vx;ux;vy;uy;vz;uz;ax;ay;az];

z_fun = [];

xdot = SX.sym('xdot',nx,1);
impl_f = xdot - x_dot;

%% Objectives and constraints

% inner objectives
h = [px(n);py(n);pz(n);vx;vy;vz;ux;uy;uz];
hN = h(1:nyN);

% outer objectives
obji = 0.5*(h-refs)'*diag(Q)*(h-refs);
objN = 0.5*(hN-refN)'*diag(QN)*(hN-refN);

obji_GGN = 0.5*(aux-refs)'*(aux-refs);
objN_GGN = 0.5*(auxN-refN)'*(auxN-refN);

% general inequality path constraints
general_con = []; 
general_con_N = []; 

%% NMPC discretizing time length [s]

Ts_st = 0.2; % shooting interval time

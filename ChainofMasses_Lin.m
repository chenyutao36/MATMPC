%------------------------------------------%
% Chain of masses connected by linear springs

% from "Efficient direct multiple shooting for nonlinear model predictive
% control on long horizons", Kirches, 2012

%------------------------------------------%

%% Dimensions

n=15;
nx=n*3+(n-1)*3;
nu=3;
np=0;
ny=3*(n-1)+3+3;
nyN=3*(n-1)+3;
nc=3;
ncN=nx;

p0x=0;p0y=0;p0z=0;

import casadi.*

states   = SX.sym('states',nx,1);
controls = SX.sym('controls',nu,1);
params   = SX.sym('paras',np,1);
refs     = SX.sym('refs',ny,1);
refN     = SX.sym('refs',nyN,1);
Q        = SX.sym('Q',ny,ny);
QN       = SX.sym('QN',nyN,nyN);

k=0.1;
lr=0.55;
m=0.45;
g=9.81;

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

%% Objectives and constraints

h = [px(n);py(n);pz(n);vx;vy;vz;ux;uy;uz];
% h = [theta;du];

hN = h(1:nyN);

h_fun=Function('h_fun', {states,controls,params}, {h},{'states','controls','params'},{'h'});
hN_fun=Function('hN_fun', {states,params}, {hN},{'states','params'},{'hN'});

ineq=[ux;uy;uz];
ineqN=states;

ineq_fun=Function('ineq_fun', {states,controls,params}, {ineq},{'states','controls','params'},{'ineq'});
ineqN_fun=Function('ineqN_fun', {states,params}, {ineqN},{'states','params'},{'ineqN'});

lb_ineq=SX.sym('lb_ineq',length(ineq),1);
ub_ineq=SX.sym('ub_ineq',length(ineq),1);
lbN_ineq=SX.sym('lbN_ineq',length(ineqN),1);
ubN_ineq=SX.sym('ubN_ineq',length(ineqN),1);
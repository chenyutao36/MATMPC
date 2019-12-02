%------------------------------------------%

% Two-stage turbocharged gasoline engine airpath

% T. Albin, D. Ritter, D. Abel, N. Liberda, R. Quirynen and M. Diehl, "Nonlinear MPC for a two-stage turbocharged gasoline engine airpath",
% 2015 54th IEEE Conference on Decision and Control (CDC), Osaka, 2015, pp. 849-856.
% doi: 10.1109/CDC.2015.7402335


% typical configuration: N=30, Ts=Ts_st=0.05, shifting, RTI=no

%------------------------------------------%


%% Dimensions

nx = 4;  % No. of differential states
nu = 2;  % No. of controls
nz = 2;  % No. of algebraic states
ny = 3; % No. of outputs
nyN = 1; % No. of outputs at the terminal point
np = 2; % No. of model parameters
nc = 3; % No. of general constraints
ncN = 3; % No. of general constraints at the terminal point
nbx = 2; % No. of bounds on states
nbu = 2; % No. of bounds on controls

% state and control bounds
nbx_idx = 3:4; % indexs of states which are bounded
nbu_idx = 1:2; % indexs of controls which are bounded

%% create variables

import casadi.*

states   = SX.sym('states',nx,1);   % differential states
controls = SX.sym('controls',nu,1); % control input
alg      = SX.sym('alg',nz,1);      % algebraic states
params   = SX.sym('paras',np,1);    % parameters
refs     = SX.sym('refs',ny,1);     % references of the first N stages
refN     = SX.sym('refs',nyN,1);    % reference of the last stage
Q        = SX.sym('Q',ny,1);        % weighting matrix of the first N stages
QN       = SX.sym('QN',nyN,1);      % weighting matrix of the last stage
aux      = SX.sym('aux',ny,1);      % auxilary variable
auxN     = SX.sym('auxN',nyN,1);    % auxilary variable

%% Dynamics
c1 = 25;
c2 = 0.004;
c3 = 8e3;
c4 = 0.5;
c5 = 40;
c6 = 1e-2;
c7 = 4e3;
c8 = 1;

ahp = 12215930.9124;
bhp = -9211942.0713;
alp = 5165385.2243;
blp = -4728474.4906;

a = [0, 1];
b = [1, -1];
g = [65, 5];
p = [1.5, 0.04];

x1 = states(1);
x2 = states(2);
u1 = states(3);
u2 = states(4);
du1 = controls(1);
du2 = controls(2);
z1 = alg(1);
z2 = alg(2);
nICE = params(1);        % engine speed
dist = params(2);        % disturbance

z1Smooth = smoothen(z1);
z2Smooth = smoothen(z2);

u1Smooth = smooth_fun(x1*x2, p, a)*smooth_fun(u1, g, b);
u2Smooth = 1-(u2/100);

% explicit ODE RHS
x_dot=[c1*(z1^(1.5) - z1^(1.25))*sqrt(smoothen(z1Smooth^(-1.5) - z1Smooth^(-1.75))) - c2*nICE*x2*(x1^(1.29) - x1);
       c5*z1*(z2^(1.5) - z2^(1.25))*sqrt(smoothen(z2Smooth^(-1.5) - z2Smooth^(-1.75))) - c6*nICE*x1*(x2^(1.29) - x2);
       du1;
       du2];

% algebraic function
z_fun = [-x1*x2 + c3/nICE*sqrt(smoothen(z1Smooth^(0.5) - z1Smooth^(0.25)))*(z1Smooth^(0.5) + c4*u1Smooth);
         -x1*x2 + c7/nICE*z1*sqrt(smoothen(z2Smooth^(0.5) - z2Smooth^(0.25)))*(z2Smooth^(0.5) + c8*u2Smooth)];

% implicit ODE: impl_f = 0
xdot = SX.sym('xdot', nx,1);
impl_f = xdot - x_dot;
     
%% Objectives and constraints

% inner objectives
h = [x1 * x2 + dist;
     du1;
     du2];

hN = h(1:nyN);

% outer objectives
obji = 0.5*(h-refs)'*diag(Q)*(h-refs);
objN = 0.5*(hN-refN)'*diag(QN)*(hN-refN);

obji_GGN = 0.5*(aux-refs)'*(aux-refs);
objN_GGN = 0.5*(auxN-refN)'*(auxN-refN);

% general inequality constraints
general_con = [x1 * x2 + dist;
               smoothen(sqrt(alp * x1 + blp));
               smoothen(sqrt(ahp * x2 + bhp))];
general_con_N = [x1 * x2 + dist;
                 smoothen(sqrt(alp * x1 + blp));
                 smoothen(sqrt(ahp * x2 + bhp))];

%% NMPC discretizing time length [s]

Ts_st = 0.05; % shooting interval time

%% Auxilliary Functions
function [ out ] = smooth_fun( x, p, a )
out = a(1) + a(2) ./ (1 + exp(-(x-p(1))/p(2)));
end

function [ out ] = smoothen( in )
EPS = 1e-4;
out = 1/2*sqrt(in.^2 + EPS) + 1/2.*in;
end

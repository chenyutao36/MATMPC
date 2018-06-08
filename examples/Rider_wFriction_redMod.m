% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File Name: casadi_mex_generation.m
%
% Designed By: Andrea De Simoi, Yutao Chen (DEI - unipd)
% Company    : VI-grade
% Project    : VI-Rider MPC fo HRC
% Purpose    : Code generation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Dimensions

srb_model_states_enum_redMod;
srb_model_controls_enum;
srb_model_parameters_enum;
mpc_controller_sizes;

%sizes definition
nx = evStatesSize; % number of states: see related enum
nu = evCtlSize;    % number of controls: see related enum
ny = cfLength;     % number of elements of che CF
nyN= cfLengthN;    % number of elements of che CF in the last point of the preview
np = evParSize;    % number of parameters: see related enum
nc = 0; % No. of general constraints
ncN= 0; % No. of general constraints at the terminal point
nbx= 2; % No. of bounds on states
nbu= 0; % No. of bounds on controls

% state and control bounds
nbx_idx = [evThr;evBrk]; % indexs of states which are bounded
nbu_idx = 0; % indexs of controls which are bounded

%% create variables

import casadi.*

states   = SX.sym('states',nx,1);
controls = SX.sym('controls',nu,1);
params   = SX.sym('paras',np,1);
refs     = SX.sym('refs',ny,1);
refN     = SX.sym('refs',nyN,1);
Q        = SX.sym('Q',ny,ny);
QN       = SX.sym('QN',nyN,nyN);

m      = params(evM);                          % vehicle total mass
g      = params(evG);                          % gravitational acceleration
p      = params(evP);                          % vehicle wheel base
b      = params(evB);                          % vehicle rear to cm distance    
h      = params(evHcm);                        % vehicle total CM height wrt ground
Izz    = params(evIzz);                        % vehicle local CM inertia zz 
Iyy    = params(evIyy);                        % vehicle local CM inertia yy
Ixx    = params(evIxx);                        % vehicle local CM inertia xx
Rwf    = params(evRwf);                        % front wheel radius
Rwr    = params(evRwr);                        % rear  wheel radius   
caster = params(evCaster);                     % front fork caster angle                     
fw_brake_trq_max  = params(evFwBrkTrqMax);     % max brake torque front wheel
rw_brake_trq_max  = params(evRwBrkTrqMin);     % max brake torque rear wheel
fw_brake_trq_bias = params(evFwBrkTrqBias);    % front/rear brake ripartition (bias) 
air_density       = params(evAirDensity);      % air density
drag_coeff_acc    = params(evDragCoeffAcc);    % drag coefficient
front_section_area= params(evFrontSecArea);    % vehicle front section area  
% tire_coeffs_f1    = params(18);               % front tire lateral slip gain coefficient [parameter used by linear tire version]
% tire_coeffs_f2    = params(19);               % front tire roll angle gain coefficient   [parameter used by linear tire version]
% tire_coeffs_r1    = params(20);               % rear tire lateral slip gain coefficient  [parameter used by linear tire version] 
% tire_coeffs_r2    = params(21);               % rear tire roll angle gain coefficient    [parameter used by linear tire version]
engTM      = params(evEngTrqMax);              % engine instantaneous max torque
engTm      = params(evEngTrqMin);              % engine instantaneous min torque
gratio     = params(evGearRatio);              % transmission instantaneous gear ratio
prm        = params(evPrm);                    % transmission primary ratio  
fin        = params(evFin);                    % transmission final ratio
LTscaling  = params(evLTScaling);              % Load Transfer scaling factor
% Fyfscaling = params(28);                      % Front tire Y force scaling factor  [parameter used by linear tire version] 
% Fyrscaling = params(29);                      % Rear  tire Y force scaling factor	 [parameter used by linear tire version]
% psi_s      = params(30);                      % Instantaneous static yaw           [debug parameter]
psi_s_first       = params(evPsiSFirst);       % Instantaneous path curvature
theta_lim         = params(evTirLimTheta);     % Limit inclination angle due to tires
barrier_order     = params(evIncBarrOrder);    % Order of thera barrier (only even usable).
Fzf_fadeout_k     = params(evFzfFoShapeFactor);% front Fz fade out shape factor 
Fzf_fadeout_theta = params(evFzfFoLimTheta);   % front Fz fade out inclination angle limit (at which fo gain = 0.5)
Fzr_fadeout_k     = params(evFzrFoShapeFactor);% rear  Fz fade out shape factor
Fzr_fadeout_theta = params(evFzrFoLimTheta);   % rear  Fz fade out inclination angle limit (at which fo gain = 0.5)
Fyf_fadeout_k     = params(evFyfFoShapeFactor);% front Fy fade out shape factor
Fyf_fadeout_lambda= params(evFyfFoLimLambda);  % front Fy fade out lateral slip limit (at which fo gain = 0.5)
Fyr_fadeout_k     = params(evFyrFoShapeFactor);% rear  Fy fade out shape factor
Fyr_fadeout_lambda= params(evFyrFoLimLambda);  % rear  Fy fade out lateral slip limit (at which fo gain = 0.5)
Fzf_m = params(evMinFzRefF);                   % reference normal force used to generate min Fyf polynomial
Fzf_M = params(evMaxFzRefF);                   % reference normal force used to generate max Fyf polynomial
Fzr_m = params(evMinFzRefR);                   % reference normal force used to generate min Fyr polynomial
Fzr_M = params(evMaxFzRefR);                   % reference normal force used to generate max Fyr polynomial

mu    = 0.9;                                   % friction coeff
theta_lim = theta_lim*mu;                      % adapt the max roll 

srb_tire_poly_constants;

e_psi = states(evEpsi);   % Yaw tracking error wrt reference trajectory  
e_y   = states(evEy);     % Path distance wrt reference trajectory
% s     = states(evS);      % Archlength instantaneous position
% t     = states(evT);      % Time
% x     = states(evX);      % Gyro global x position ISO reference frame  
% y     = states(evY);      % Gyro global y position ISO reference frame 
yaw   = states(evYaw);    % Gyro global yaw 
theta = states(evTheta);  % CM roll angle wrt GYRO
vx    = states(evVx);     % Gyro global vx  
vy    = states(evVy);     % Gyro global vy  
yr    = states(evYr);     % Gyro global yaw rate 
rr    = states(evRr);     % CM roll rate wrt GYRO
delta = states(evDelta);  % steering angle 
thr   = states(evThr);    % throtle        [0-1]
brk   = states(evBrk);    % brake          [0-1]
pitch = states(evPitch);  % Sentinel variable (simil-pitch)
pr    = states(evPr);     % Sentinel variable dot (simil-pitch rate) 

delta_dot = controls(evDeltaDot); % Steering derivative    
thr_dot   = controls(evThrDot);   % Throttle derivative          
brk_dot   = controls(evBrkDot);   % Brrake   derivetive 

%% Intermediate States
delta_G  = atan(tan(delta)*cos(caster)/...
              (cos(theta)+tan(delta)*sin(theta)*sin(caster)));

% Inertia values correction due to solver instability. 
% TODO: Investigate why they are required. Missing full inertia tensor?
Ixx_cor  = Ixx*20;
Izz_cor  = Izz*20;

% tire profile approximatet with 0.06 front and 0.09 rear
st = -1*sin(theta);
ct = 1-cos(theta);
% Rear Contact Patch position 
CPr_x    = 0.0;
CPr_y    = 0.09*st ; 
CPr_z    = 0.09*ct; 
% Front Contact Patch position
CPf_x    = p + 0.06*st*sin(delta_G);
CPf_y    = 0.06*st*cos(delta_G);
CPf_z    = 0.06*ct   ;
% CM position according to Contact patches positioning
CM_y     = -h*sin(theta);
CM_z     = h*cos(theta);
CM_dy_gr = (b*CPf_y + (CPf_x-b)*CPr_y)/CPf_x ;  
CM_dz_gr = (b*CPf_z + (CPf_x-b)*CPr_z)/CPf_x ;
% update hcm and theta according to CM position and contact patch positions
hcm        = sqrt((CM_y-CM_dy_gr)^2 + power((CM_z+CM_dz_gr),2));
thetaIdeal = -1*atan((CM_y-CM_dy_gr)/(CM_z+CM_dz_gr));
sideSlipPf = atan((CPf_y-CPr_y)/(CPf_x-CPr_x));

% vx and vy are already in gyro reference frame 
vx_gyro  = vx ;  %( cos(yaw)*vx + sin(yaw)*vy);
vy_gyro  = vy ;  %(-sin(yaw)*vx + cos(yaw)*vy);
% correct vx gyro in function of the side slip angle caused by tire profiles
% vx_gyro  = ( cos(sideSlipPf)*vx - sin(sideSlipPf)*vy );
% vy_gyro  = ( sin(sideSlipPf)*vx + cos(sideSlipPf)*vy );
v_module = sqrt(vx*vx + vy*vy) ; 

s_dot = (1.0/(1.0-psi_s_first*e_y))*(vx_gyro*cos(e_psi+sideSlipPf) - vy_gyro*sin(e_psi+sideSlipPf)) ;

% front/rear wheel lateral slip
lambda_f = atan((vy_gyro + p*yr)/vx_gyro) - delta_G ;
lambda_r = atan(vy_gyro/vx_gyro);

% transmission torque min/max calculation
rw_trans_trq_max  = engTM/gratio/prm/fin ;
rw_trans_trq_min  = engTm/gratio/prm/fin ;
% drive 
driveTorque = thr*(rw_trans_trq_max - rw_trans_trq_min)+rw_trans_trq_min;

% %%%%%%%%%%%%%%%%%%%%%%% %
%   Fz Fx Fy calculation  %
% %%%%%%%%%%%%%%%%%%%%%%% %

% static Fx first implementation
% front/rear tire X force in tire ISO reference frame
% Fxf = -brk*fw_brake_trq_bias*fw_brake_trq_max/Rwf;
% Fxr = ( driveTorque - brk*(1.0-fw_brake_trq_bias)*rw_brake_trq_max)/Rwr;

% Fx calculation with saturation accordinf to the sentinel variable
Fxf =  (-brk*fw_brake_trq_bias*fw_brake_trq_max/Rwf)*-1*(Fz_sat(-pitch/0.05)-1);
Fxr =  (( driveTorque - brk*(1.0-fw_brake_trq_bias)*rw_brake_trq_max)/Rwr)*-1*(Fz_sat(pitch/0.05)-1);

% friction
Fxf = Fxf*mu;
Fxr = Fxr*mu;

% %%%%%%%%%%%%%%%%%%%%%%%%%% %
%   Load transfer estimation  %%
%  %%%%%%%%%%%%%%%%%%%%%%%%%% %%

% First Load transfer estimation impl. time/space invariant (STATIC LOAD TRANSFER)
LT = (Fxf+Fxr)*h*cos(theta)/p*LTscaling ;  %% LTscaling scale factor usable to disable load transfer

% Not normalized Fzf/r (first implementation)
% Fzf = m*g*b/p - LT;
% Fzr = m*g*(p-b)/p + LT; 

% Apply LT to static normal forces and normalize to m*g
Fzf_drop = 1/(1+exp(-Fzf_fadeout_k*(theta+Fzf_fadeout_theta)))... 
          	 *1/(1+exp( Fzf_fadeout_k*(theta-Fzf_fadeout_theta))) ; 
Fzr_drop = 1/(1+exp(-Fzr_fadeout_k*(theta+Fzr_fadeout_theta)))... 
             *1/(1+exp( Fzr_fadeout_k*(theta-Fzr_fadeout_theta))) ;
FzfNorm  = (m*g*b/p     - LT)/(m*g)*Fzf_drop ;
FzrNorm  = (m*g*(p-b)/p + LT)/(m*g)*Fzr_drop ;
% saturate normalized Fzf/r and get back the the value in the range [0 - mg]
% front/rear tire normal force in tire ISO reference frame
Fzf =  Fz_sat(FzfNorm)*(m*g)*-1*(Fz_sat(-pitch/0.05)-1) ;
Fzr =  Fz_sat(FzrNorm)*(m*g)*-1*(Fz_sat( pitch/0.05)-1) ; 

% evaluate tire surfaces
srb_tire_poly_surfaces;

Fyf = Fyf_m + (Fyf_M-Fyf_m)*(Fzf-Fzf_m)/(Fzf_M-Fzf_m) ;
Fyr = Fyr_m + (Fyr_M-Fyr_m)*(Fzr-Fzr_m)/(Fzr_M-Fzr_m) ;

satslip_f = -1*(Fz_sat(-pitch/0.05)-1)*1/(1+exp(-Fyf_fadeout_k*(lambda_f+Fyf_fadeout_lambda)))... 
                *1/(1+exp( Fyf_fadeout_k*(lambda_f-Fyf_fadeout_lambda))) ;
satslip_r = -1*(Fz_sat( pitch/0.05)-1)*1/(1+exp(-Fyr_fadeout_k*(lambda_r+Fyr_fadeout_lambda)))... 
                *1/(1+exp( Fyr_fadeout_k*(lambda_r-Fyr_fadeout_lambda))) ;
Fyf = Fyf*satslip_f ;
Fyr = Fyr*satslip_r ;

% friction
Fyf = Fyf*mu;
Fyr = Fyr*mu;


% linear tire version
% Fyf = (tire_coeffs_f1*lambda_f + tire_coeffs_f2*theta)*Fzf;
% Fyr = (tire_coeffs_r1*lambda_r + tire_coeffs_r2*theta)*Fzf;

% drag force
Fd  = -0.5*air_density*drag_coeff_acc*front_section_area*vx_gyro*vx_gyro ;
% centrifugal force
F_c = -m*vx_gyro*yr ;

%front/rear wheel omega
omega_wf =  vx_gyro/Rwf ; 
omega_wr =  vx_gyro/Rwr ;

rollBarrier = power((theta/theta_lim),barrier_order) ;
% hardcoded sentinel limit value
pitchBarrier= power((pitch/0.02),barrier_order);

%% Model Equations
srb_c_equations_wrt_gyro_rhs;

x_dot= [ (yr/s_dot) - psi_s_first                ; ... % e_psi
(vx_gyro*sin(e_psi) + vy_gyro*cos(e_psi))/s_dot ; ... % e_y
% 1                                               ; ... % s
% 1/s_dot                                         ; ... % t
% (vx*cos(yaw)-vy*sin(yaw))/s_dot                 ; ... % x     vx,vy in global (ISO))
% (vx*sin(yaw)+vy*cos(yaw))/s_dot                 ; ... % y     ref frame
yr/s_dot                                        ; ... % yaw
rr/s_dot                                        ; ... % theta
vx_dot/s_dot                                    ; ... % vx
vy_dot/s_dot                                    ; ... % vy
yr_dot/s_dot                                    ; ... % yr
rr_dot/s_dot                                    ; ... % rr
delta_dot/s_dot                                 ; ... % delta
thr_dot/s_dot                                   ; ... % thr
brk_dot/s_dot                                   ; ... % brk
pr/s_dot                                        ; ... % pitch
pr_dot/s_dot 
];    % pr   

xdot = SX.sym('xdot',nx,1);
impl_f = xdot - x_dot;

%%intermediate states dump
intStates = [delta_G, vx_gyro, vy_gyro, lambda_f, lambda_r, Fxf, Fxr,...
             Fzf, Fzr, Fyf, Fyr, vx_dot, vy_dot, yr_dot, rr_dot, hcm,...
             thetaIdeal, sideSlipPf, pitch, pr, pr_dot];

intermStates  = Function('intermStates', {states,controls,params}, {intStates},{'states','controls','params'},{'intStates'});            

%% Objectives and constraints

h = [states(evEpsi) + atan( (states(evVy)+params(evB)*states(evYr))/states(evVx) );...
     states(evEy);...
     sqrt( states(evVx)^2+states(evVy)^2 );...
     power( (states(evTheta)/params(evTirLimTheta)), params(evIncBarrOrder) );...
     power( (states(evPitch)/0.02), params(evIncBarrOrder) );...
     states(evBrk)*states(evTheta);...
     controls(evDeltaDot);...
     controls(evThrDot);...
     controls(evBrkDot);...
     states(evThr)*states(evBrk)];

hN = h(1:nyN);

% general inequality constraints
general_con = [];
general_con_N = [];

%% NMPC sampling time [s]

Ts = 1; % simulation sample time
Ts_st = 1; % shooting interval time

%% build casadi function (don't touch)
h_fun=Function('h_fun', {states,controls,params}, {h},{'states','controls','params'},{'h'});
hN_fun=Function('hN_fun', {states,params}, {hN},{'states','params'},{'hN'});

path_con_fun=Function('path_con_fun', {states,controls,params}, {general_con},{'states','controls','params'},{'general_con'});
path_con_N_fun=Function('path_con_N_fun', {states,params}, {general_con_N},{'states','params'},{'general_con_N'});


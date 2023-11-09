clear all;clc;
disp( ' ' );
disp( 'MATMPC -- A (MAT)LAB based Model(M) Predictive(P) Control(C) Package.' );
disp( 'Copyright (C) 2016-2019 by Yutao Chen, University of Padova' );
disp( 'All rights reserved.' );
disp( ' ' );
disp( 'MATMPC is distributed under the terms of the' );
disp( 'GNU General Public License 3.0 in the hope that it will be' );
disp( 'useful, but WITHOUT ANY WARRANTY; without even the implied warranty' );
disp( 'of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.' );
disp( 'See the GNU General Public License for more details.' );
disp( ' ' );
disp( ' ' );
disp('---------------------------------------------------------------------------------');

%% Insert Model here
addpath([pwd,'/examples']);

settings.model='InvertedPendulum'; % see the folder "examples" for details

run(settings.model);

%%
% The code reported in this script refers to the papers reported below
% - [1] Y. Chen, M. Bruschetta, E. Picotti and A. Beghi, "MATMPC - A MATLAB 
% Based Toolbox for Real-time Nonlinear Model Predictive Control," 2019 18th 
% European Control Conference (ECC), 2019, pp. 3365-3370, doi: 10.23919/ECC.2019.8795788.
% - [2] Y. Chen, M. Bruschetta, D. Cuccato and A. Beghi, "An Adaptive Partial 
% Sensitivity Updating Scheme for Fast Nonlinear Model Predictive Control," in 
% IEEE Transactions on Automatic Control, vol. 64, no. 7, pp. 2712-2726, July 2019, 
% doi: 10.1109/TAC.2018.2867916.

import casadi.*

lambda=SX.sym('lambda',nx,1);            % the i th multiplier for equality constraints
mu_u=SX.sym('mu_u',nu,1);                  % the i th multiplier for bounds on controls
mu_x=SX.sym('mu_x',nbx,1);                  % the i th multiplier for bounds on controls
mu_g=SX.sym('mu_g',nc,1);                  % the i th multiplier for bounds on controls
muN_g=SX.sym('muN_g',ncN,1);                 % the N th multiplier for inequality constraints

%% Generate some functions

% Computing the Jacobian of the implicit function [1, eqs.(1a)--(1e)]. A
% "lighther" description is also reported in [2, eqs.(1a)--(1e)]. For this
% part refers to the online documentation of CaSADI: https://web.casadi.org/docs/#document-function
f_fun  = Function('f_fun', {states,controls,params,alg}, {SX.zeros(nx,1)+x_dot},{'states','controls','params','alg'},{'xdot'});

% Initializate the dimension of the symbolic vectors in CaSADI
impl_jac_f_x = SX.zeros(nx,nx)+jacobian(impl_f,states);
impl_jac_f_u = SX.zeros(nx,nu)+jacobian(impl_f,controls);
impl_jac_f_xdot = SX.zeros(nx,nx)+jacobian(impl_f,xdot);

% Allocating the Functions in CaSADI. For this part refers to the online
% documentation: https://web.casadi.org/docs/#document-function
impl_f_fun = Function('impl_f_fun',{states,controls,params,xdot,alg},{SX.zeros(nx,1) + impl_f});
impl_jac_f_x_fun = Function('impl_jac_f_x_fun',{states,controls,params,xdot,alg},{impl_jac_f_x});
impl_jac_f_u_fun = Function('impl_jac_f_u_fun',{states,controls,params,xdot,alg},{impl_jac_f_u});
impl_jac_f_xdot_fun = Function('impl_jac_f_xdot_fun',{states,controls,params,xdot,alg},{impl_jac_f_xdot});

% In case the number of algebraic states is greater than zero
if nz>0
    % F
    g_fun  = Function('g_fun', {states,controls,params,xdot,alg}, {SX.zeros(nz,1)+z_fun},{'states','controls','params','xdot','alg'},{'zfun'});
    
    % Initializate the dimension of the symbolic vectors in CaSADI in the
    % case the number of algebraic states is greater than zero
    impl_jac_f_z = SX.zeros(nx,nz)+jacobian(impl_f,alg);
    impl_jac_g_x = SX.zeros(nz,nx)+jacobian(z_fun,states);
    impl_jac_g_u = SX.zeros(nz,nu)+jacobian(z_fun,controls);
    impl_jac_g_z = SX.zeros(nz,nz)+jacobian(z_fun,alg);
    
    % Allocating the Functions in CaSADI in the case the number of
    % algebraic states is greater than zero. For this part refers to the
    % online documentation of CaSADI: https://web.casadi.org/docs/#document-function
    impl_jac_f_z_fun = Function('impl_jac_f_z_fun',{states,controls,params,xdot,alg},{impl_jac_f_z});
    impl_jac_g_x_fun = Function('impl_jac_g_x_fun',{states,controls,params,xdot,alg},{impl_jac_g_x});
    impl_jac_g_u_fun = Function('impl_jac_g_u_fun',{states,controls,params,xdot,alg},{impl_jac_g_u});
    impl_jac_g_z_fun = Function('impl_jac_g_z_fun',{states,controls,params,xdot,alg},{impl_jac_g_z});
else
    % In the case the number of algrebraic states is equal to zero. In
    % CaSADI the expression {0} refers to an empty function. For more
    % information refers to the online documentation of CaSADI:
    % https://web.casadi.org/docs/#document-function
    g_fun  = Function('g_fun', {states,controls,params,xdot,alg}, {0},{'states','controls','params','xdot','alg'},{'zfun'});
    
    % Allocating the Functions in CaSADI in the case the number of
    % algebraic states is equal to zero. For more info refers to the online
    % documentation of CaSADI: https://web.casadi.org/docs/#document-function
    impl_jac_f_z_fun = Function('impl_jac_f_z_fun',{states,controls,params,xdot,alg},{0});
    impl_jac_g_x_fun = Function('impl_jac_g_x_fun',{states,controls,params,xdot,alg},{0});
    impl_jac_g_u_fun = Function('impl_jac_g_u_fun',{states,controls,params,xdot,alg},{0});
    impl_jac_g_z_fun = Function('impl_jac_g_z_fun',{states,controls,params,xdot,alg},{0});
end

Sx = SX.sym('Sx',nx,nx);
Su = SX.sym('Su',nx,nu);
% In the case the number of algebraic states is equal to zero
if nz==0    
    % This part implements the ADJ-RTI scheme (see [2, Sec.II]. For the
    % "jacobian" function refers to
    % https://web.casadi.org/docs/#calculus-algorithmic-differentiation.
    % Instead, the "jtimes" function is described in the subsection 3.9.1.
    % of the online documentation of CaSADI. 
    % 
    % jtimes = For calculating a Jacobian-times-vector product, the jtimes function 
    % – performing forward mode AD – is often more efficient than creating the full 
    % Jacobian and performing a matrix-vector multiplication:
    vdeX = SX.zeros(nx,nx);
    vdeX = vdeX + jtimes(x_dot,states,Sx);
    vdeU = SX.zeros(nx,nu) + jacobian(x_dot,controls);
    vdeU = vdeU + jtimes(x_dot,states,Su);
    vdeFun = Function('vdeFun',{states,controls,params,Sx,Su,alg},{vdeX,vdeU});
    
    % ERK stands for explicit Runge–Kutta integrator (see [1, Sec. V])
    adjW = SX.zeros(nx+nu,1) + jtimes(x_dot, [states;controls], lambda, true);
    adj_ERK_fun = Function('adj_ERK_fun',{states,controls,params,lambda,alg},{adjW});
% in the opposite case (i.e., the number of algrebraic states is greater
% than zero)
else
    % ERK stands for explicit Runge–Kutta integrator (see [1, Sec. V])
    vdeFun = Function('vdeFun',{states,controls,params,Sx,Su,alg},{0,0});
    adj_ERK_fun = Function('adj_ERK_fun',{states,controls,params,lambda,alg},{0});
end
    
%% objective and constraints

% The inner function and the path constraint function are described in [2].
% For more information regarding the meaning, please refer to the following
% papers:
% - [3] C.Kirches, L. Wirsching, H. G. Bock, J. P. Schlöder, "Efficient direct 
% multiple shooting for nonlinear model predictive control on long
% horizons", Journal of Process Control, Volume 22, Issue 3, March 2012,
% Pages 540--550, doi: 10.1016/j.jprocont.2012.01.008
% - [4] R. Quirynen, M. Vukov, M. Zanon, M. Diehl, "Autogenerating microsecond 
% solvers for nonlinear MPC: A tutorial using ACADO integrators", Optimal
% Control, doi: https://doi.org/10.1002/oca.2152
% - [5] E. Rossi, M. Bruschetta, R. Carli, Y. Chen and M. Farina, "Online Nonlinear 
% Model Predictive Control for tethered UAVs to perform a safe and constrained maneuver,"
% 2019 18th European Control Conference (ECC), 2019, pp. 3996-4001, doi: 10.23919/ECC.2019.8796032.
% - [6] T. Albin, D. Ritter, D. Abel, N. Liberda, R. Quirynen and M. Diehl, 
% "Nonlinear MPC for a two-stage turbocharged gasoline engine airpath", 2015 
% 54th IEEE Conference on Decision and Control (CDC), Osaka, 2015, pp. 849-856.
% doi: 10.1109/CDC.2015.7402335
% Some numerical examples are also presented in [1]
h_fun=Function('h_fun', {states,controls,params}, {h},{'states','controls','params'},{'h'});
hN_fun=Function('hN_fun', {states,params}, {hN},{'states','params'},{'hN'});
path_con_fun=Function('path_con_fun', {states,controls,params}, {general_con},{'states','controls','params'},{'general_con'});
path_con_N_fun=Function('path_con_N_fun', {states,params}, {general_con_N},{'states','params'},{'general_con_N'});

% This part refers to [2, eq.(4)]. 
gxi = jacobian(obji,states)' + SX.zeros(nx,1);
gui = jacobian(obji,controls)' + SX.zeros(nu,1);
gxN = jacobian(objN,states)' + SX.zeros(nx,1);

% Gauss-Newton Hessian approximation [2, eq.(5)]. 
Hz = hessian(obji_GGN, aux) + SX.zeros(ny, ny);
HzN = hessian(objN_GGN, auxN) + SX.zeros(nyN, nyN);

obj_vec = sqrt(diag(Q))*h;
objN_vec = sqrt(diag(QN))*hN;
Jxi = jacobian(obj_vec, states) + SX.zeros(ny, nx);
Jui = jacobian(obj_vec, controls) + SX.zeros(ny, nu);
JxN = jacobian(objN_vec, states) + SX.zeros(nyN, nx);

Cxi = jacobian(general_con, states) + SX.zeros(nc, nx);
Cui = jacobian(general_con, controls) + SX.zeros(nc, nu);
CxN = jacobian(general_con_N, states) + SX.zeros(ncN, nx);

% Objective functions at the terminal stage (objN_fun) and not)
obji_fun = Function('obji_fun',{states,controls,params,refs,Q},{obji+SX.zeros(1,1)});
objN_fun = Function('objN_fun',{states,params,refN,QN},{objN+SX.zeros(1,1)});

% Creating functions in CaSADI
Hz_fun = Function('Hz_fun',{aux,params,refs,Q},{Hz});
HzN_fun = Function('HzN_fun',{auxN,params,refN,QN},{HzN});
Hw = Hz_fun(h,params,refs,Q) + SX.zeros(ny,ny);
HwN = HzN_fun(hN,params,refN,QN)+ SX.zeros(nyN,nyN);

% Creating functions in CaSADI
Hi_fun=Function('Hi_fun',{states,controls,params,refs,Q},{Hw});
HN_fun=Function('HN_fun',{states,params,refN,QN},{HwN});

Ji_fun=Function('Ji_fun',{states,controls,params,refs,Q},{Jxi,Jui});
JN_fun=Function('JN_fun',{states,params,refN,QN},{JxN});

gi_fun=Function('gi_fun',{states,controls,params,refs,Q},{gxi, gui});
gN_fun=Function('gN_fun',{states,params,refN,QN},{gxN});

Ci_fun=Function('Ci_fun',{states,controls,params},{Cxi, Cui});
CN_fun=Function('CN_fun',{states,params},{CxN});


dobj = [gxi;gui];
dobjN = gxN;
    
Jxb = zeros(nbx,nx+nu);
for i=1:nbx
    Jxb(i,nbx_idx(i)) = 1.0; 
end

Jub = zeros(nu,nx+nu);
for i=1:nu
    Jub(i,nx+i)=1.0;
end
adj_dB = SX.zeros(nx+nu,1) + Jxb'*mu_x + Jub'*mu_u;
if nc>0
   adj_dB = adj_dB + jtimes(general_con, [states;controls], mu_g, true);
end

JxbN = zeros(nbx,nx);
for i=1:nbx
    JxbN(i,nbx_idx(i)) = 1.0; 
end
adj_dBN = SX.zeros(nx,1) + JxbN'*mu_x;
if ncN>0
    adj_dBN = adj_dBN + jtimes(general_con_N, states, muN_g, true);
end

adj_fun = Function('adj_fun',{states,controls,params,refs,Q,lambda,mu_x,mu_u,mu_g},{dobj, adj_dB});
adjN_fun = Function('adjN_fun',{states,params,refN, QN, mu_x,muN_g},{dobjN, adj_dBN});

%% Code generation and Compile

generate=input('Would you like to generate the source code?(y/n)','s');

if strcmp(generate,'y')

    disp('                           ');
    disp('    Generating source code...');
    
    cd model_src
      
    opts = struct( 'main', false, 'mex' , true ) ; 
    h_fun.generate('h_fun.c',opts);
    f_fun.generate('f_fun.c',opts);
    g_fun.generate('g_fun.c',opts);
    path_con_fun.generate('path_con_fun.c',opts);
    path_con_N_fun.generate('path_con_N_fun.c',opts);
   
    opts = struct('main',false,'mex',false,'with_header',true);
    cd ../mex_core
        P = CodeGenerator ('casadi_src.c', opts) ;
        P.add(f_fun);
        P.add(g_fun);
        P.add(vdeFun);
        P.add(adj_ERK_fun);
        P.add(adj_fun);
        P.add(adjN_fun);       
        P.add(impl_f_fun);
        P.add(impl_jac_f_x_fun);
        P.add(impl_jac_f_u_fun);
        P.add(impl_jac_f_xdot_fun);
        P.add(impl_jac_f_z_fun);
        P.add(impl_jac_g_x_fun);
        P.add(impl_jac_g_u_fun);
        P.add(impl_jac_g_z_fun);
        P.add(h_fun);
        P.add(path_con_fun);
        P.add(path_con_N_fun);
        P.add(obji_fun);
        P.add(objN_fun);
        P.add(gi_fun);
        P.add(gN_fun);
        P.add(Hi_fun);
        P.add(HN_fun);
        P.add(Ji_fun);
        P.add(JN_fun);
        P.add(Ci_fun);
        P.add(CN_fun);
              
        P.generate();
    cd ..

disp('    Code generation completed!');

end

disp('                           ');
compile=input('Would you like to compile the source code?(y/n)','s');
if strcmp(compile,'y')
    
    disp('    Compiling...');
    
    OS_MAC = 0;
    OS_LINUX = 0;
    OS_WIN = 0;

    if ismac
        OS_MAC = 1;
    elseif isunix
        OS_LINUX = 1;
    elseif ispc
        OS_WIN = 1;
    else
        disp('    Platform not supported')
    end
    
    options = '-largeArrayDims';

    if OS_WIN
       CC_FLAGS='CXXFLAGS="$CXXFLAGS -Wall"'; % use MinGW not VS studio
    end
    if OS_LINUX 
	    CC_FLAGS = 'GCC="/usr/bin/gcc"';
    end
    if OS_MAC
       CC_FLAGS = '';
    end
    
    OP_FLAGS='-O';
    PRINT_FLAGS='-silent';
    
    cd model_src
    mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, 'path_con_fun.c');
    mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, 'path_con_N_fun.c');
    mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, 'h_fun.c');
    mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, 'f_fun.c');
    mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, 'g_fun.c');

    cd ../mex_core
    Compile_Mex;
    cd ..

disp('    Compilation completed!');

end

%% NMPC preparation

disp('                           ');
disp('Preparing the NMPC solver...');

settings.Ts_st = Ts_st;
settings.nx = nx; 
settings.nu = nu;    
settings.nz = nz;
settings.ny = ny;    
settings.nyN= nyN;    
settings.np = np;   
settings.nc = nc;
settings.ncN = ncN;
settings.nbx = nbx;
settings.nbu = nbu;
settings.nbx_idx = nbx_idx;
settings.nbu_idx = nbu_idx;

cd data
save('settings','settings');
cd ..

clear all;

disp('NMPC solver prepared! Enjoy solving...');
disp('                           ');

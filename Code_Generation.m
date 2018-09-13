clear all;clc;
disp( ' ' );
disp( 'MATMPC -- A (MAT)LAB based Model(M) Predictive(P) Control(C) Package.' );
disp( 'Copyright (C) 2016-2018 by Yutao Chen, University of Padova' );
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
import casadi.*

lambda=SX.sym('lambdai',nx,1);            % the i th multiplier for equality constraints
mu_u=SX.sym('mu_u',nu,1);                  % the i th multiplier for bounds on controls
mu_x=SX.sym('mu_x',nbx,1);                  % the i th multiplier for bounds on controls
mu_g=SX.sym('mu_g',nc,1);                  % the i th multiplier for bounds on controls
muN_g=SX.sym('muN_g',ncN,1);                 % the N th multiplier for inequality constraints

%% Explicit Runge-Kutta 4 Integrator for simulation
s  = 1; % No. of integration steps per sample interval
DT = Ts/s;
f  = Function('f', {states,controls,params}, {x_dot},{'states','controls','params'},{'xdot'});
X=states;
U=controls; 
P=params;
for j=1:s
       [k1] = f(X, U, P);
       [k2] = f(X + DT/2 * k1, U, P);
       [k3] = f(X + DT/2 * k2, U, P);
       [k4] = f(X + DT * k3, U, P);
       X=X+DT/6*(k1 +2*k2 +2*k3 +k4);
end
Simulate_system = Function('Simulate_system', {states,controls,params}, {X}, {'states','controls','params'}, {'xf'});

%% Integrator for multiple shooting
s  = 1; % No. of integration steps per shooting interval
DT = Ts_st/s;
f_fun  = Function('f_fun', {states,controls,params}, {SX.zeros(nx,1)+x_dot},{'states','controls','params'},{'xdot'});

impl_jac_x = SX.zeros(nx,nx)+jacobian(impl_f,states);
impl_jac_u = SX.zeros(nx,nu)+jacobian(impl_f,controls);
impl_jac_xdot = SX.zeros(nx,nx)+jacobian(impl_f,xdot);
impl_f = SX.zeros(nx,1) + impl_f;
impl_f_fun = Function('impl_f_fun',{states,controls,params,xdot},{impl_f});
impl_jac_x_fun = Function('impl_jac_x_fun',{states,controls,params,xdot},{impl_jac_x});
impl_jac_u_fun = Function('impl_jac_u_fun',{states,controls,params,xdot},{impl_jac_u});
impl_jac_xdot_fun = Function('impl_jac_xdot_fun',{states,controls,params,xdot},{impl_jac_xdot});

Sx = SX.sym('Sx',nx,nx);
Su = SX.sym('Su',nx,nu);
vdeX = SX.zeros(nx,nx);
vdeX = vdeX + jtimes(x_dot,states,Sx);
vdeU = SX.zeros(nx,nu) + jacobian(x_dot,controls);
vdeU = vdeU + jtimes(x_dot,states,Su);
vdeFun = Function('vdeFun',{states,controls,params,Sx,Su},{vdeX,vdeU});

X=states;
U=controls; 
P=params;
for j=1:s
       [k1] = f_fun(X, U, P);
       [k2] = f_fun(X + DT/2 * k1, U, P);
       [k3] = f_fun(X + DT/2 * k2, U, P);
       [k4] = f_fun(X + DT * k3, U, P);
       X=X+DT/6*(k1 +2*k2 +2*k3 +k4);
end
F = Function('F', {states,controls,params}, {X + SX.zeros(nx,1)});
A = jacobian(X,states) + SX.zeros(nx,nx);
B = jacobian(X,controls) + SX.zeros(nx,nu);
D = Function('D', {states,controls,params}, {A, B});

%% objective and constraints
refs     = SX.sym('refs',ny,1);     % references of the first N stages
refN     = SX.sym('refs',nyN,1);    % reference of the last stage
Q        = SX.sym('Q',ny,1);        % weighting matrix of the first N stages
QN       = SX.sym('QN',nyN,1);      % weighting matrix of the last stage

obji_vec=sqrt(diag(Q))*(h_fun(states,controls,params)-refs);
objN_vec=sqrt(diag(QN))*(hN_fun(states,params)-refN);
Jxi = jacobian(obji_vec, states) + SX.zeros(ny, nx);
Jui = jacobian(obji_vec, controls) + SX.zeros(ny, nu);
JxN = jacobian(objN_vec, states) + SX.zeros(nyN, nx);

obji = 0.5*(h_fun(states,controls,params)-refs)'*diag(Q)*(h_fun(states,controls,params)-refs);
objN = 0.5*(hN_fun(states,params)-refN)'*diag(QN)*(hN_fun(states,params)-refN);
gxi = jacobian(obji,states)' + SX.zeros(nx,1);
gui = jacobian(obji,controls)' + SX.zeros(nu,1);
gxN = jacobian(objN,states)' + SX.zeros(nx,1);

Cxi = jacobian(general_con, states) + SX.zeros(nc, nx);
Cui = jacobian(general_con, controls) + SX.zeros(nc, nu);
CxN = jacobian(general_con_N, states) + SX.zeros(ncN, nx);

obji_fun = Function('obji_fun',{states,controls,params,refs,Q},{obji+SX.zeros(1,1)});
objN_fun = Function('objN_fun',{states,params,refN,QN},{objN+SX.zeros(1,1)});

Ji_fun=Function('Ji_fun',{states,controls,params,refs,Q},{Jxi,Jui});
JN_fun=Function('JN_fun',{states,params,refN,QN},{JxN});

gi_fun=Function('gi_fun',{states,controls,params,refs,Q},{gxi, gui});
gN_fun=Function('gN_fun',{states,params,refN,QN},{gxN});

Ci_fun=Function('Ci_fun',{states,controls,params},{Cxi, Cui});
CN_fun=Function('CN_fun',{states,params},{CxN});

dobj = [gxi;gui];
dobjN = gxN;

adj_dG = SX.zeros(nx+nu,1) + jtimes(X, [states;controls], lambda, true);

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

adj_fun = Function('adj_fun',{states,controls,params,refs,Q,lambda,mu_x,mu_u,mu_g},{dobj, adj_dG, adj_dB});
adjN_fun = Function('adjN_fun',{states,params,refN, QN, mu_x,muN_g},{dobjN, adj_dBN});

adj_dG_fun = Function('adj_dG_fun',{states,controls,params,refs,Q,lambda},{dobj, adj_dG});

%% Code generation and Compile

generate=input('Would you like to generate the source code?(y/n)','s');

if strcmp(generate,'y')

    display('                           ');
    display('    Generating source code...');
    
    cd model_src
      
    opts = struct( 'main', false, 'mex' , true ) ; 
    Simulate_system.generate('Simulate_system.c',opts);
%     D.generate('D.c',opts);
    h_fun.generate('h_fun.c',opts);
    path_con_fun.generate('path_con_fun.c',opts);
    path_con_N_fun.generate('path_con_N_fun.c',opts);
    Ji_fun.generate('Ji_fun.c',opts);
    JN_fun.generate('JN_fun.c',opts);
%     intermStates.generate('intermStates',opts);
%     costFun.generate('costFun',opts);
   
    opts = struct('main',false,'mex',false,'with_header',true);
    cd ../mex_core
        P = CodeGenerator ('casadi_src.c', opts) ;
        P.add(f_fun);
        P.add(vdeFun);
        P.add(impl_f_fun);
        P.add(impl_jac_x_fun);
        P.add(impl_jac_u_fun);
        P.add(impl_jac_xdot_fun);
        P.add(F);
        P.add(D);
        P.add(h_fun);
        P.add(path_con_fun);
        P.add(path_con_N_fun);
        P.add(obji_fun);
        P.add(objN_fun);
        P.add(gi_fun);
        P.add(gN_fun);
        P.add(Ji_fun);
        P.add(JN_fun);
        P.add(Ci_fun);
        P.add(CN_fun);
        P.add(adj_fun);
        P.add(adjN_fun);
        
        P.add(adj_dG_fun);
        
        P.generate();
    cd ..

display('    Code generation completed!');

end

display('                           ');
compile=input('Would you like to compile the source code?(y/n)','s');
if strcmp(compile,'y')
    
    display('    Compiling...');
    
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
    
    OP_FLAGS='-O';
    PRINT_FLAGS='-silent';
    
    cd model_src
    mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, 'path_con_fun.c');
%     mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, 'D.c');
    mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, 'path_con_N_fun.c');
    mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, 'h_fun.c');
    mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, 'Simulate_system.c');
    mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, 'Ji_fun.c');
    mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, 'JN_fun.c');
%     mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, 'intermStates.c');
%     mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, 'costFun.c');
       
    cd ../mex_core
    Compile_Mex;
    cd ..

display('    Compilation completed!');

end

%% NMPC preparation

display('                           ');
display('Preparing the NMPC solver...');

settings.Ts = Ts;
settings.Ts_st = Ts_st;
settings.s =s ;
settings.nx = nx; 
settings.nu = nu;    
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

display('NMPC solver prepared! Enjoy solving...');
display('                           ');

clear all;clc;
disp('---------------------------------------------');
disp('MATMPC is developed by Yutao Chen, DEI, UniPD');
disp('---------------------------------------------');

%% Insert Model here
addpath([pwd,'/examples']);

settings.model='TethUAV_param'; % see the folder "examples" for details

run(settings.model);

%%
import casadi.*

lambdai=SX.sym('lambdai',nx,1);            % the i th multiplier for equality constraints
mui=SX.sym('mui',nc,1);                  % the i th multiplier for inequality constraints
muN=SX.sym('muN',ncN,1);                 % the N th multiplier for inequality constraints

%% Explicit Runge-Kutta 4 Integrator for simulation
s  = 2; % No. of integration steps per sample interval
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
s  = 2; % No. of integration steps per shooting interval
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
       [k1] = f(X, U, P);
       [k2] = f(X + DT/2 * k1, U, P);
       [k3] = f(X + DT/2 * k2, U, P);
       [k4] = f(X + DT * k3, U, P);
       X=X+DT/6*(k1 +2*k2 +2*k3 +k4);
end
% z = [states;controls];
F = Function('F', {states,controls,params}, {X + SX.zeros(nx,1)});
A = jacobian(X,states) + SX.zeros(nx,nx);
B = jacobian(X,controls) + SX.zeros(nx,nu);
D = Function('D', {states,controls,params}, {A, B});

%% objective and constraints

obji_vec=sqrt(Q)*(h_fun(states,controls,params)-refs);
objN_vec=sqrt(QN)*(hN_fun(states,params)-refN);
Jxi = jacobian(obji_vec, states) + SX.zeros(ny, nx);
Jui = jacobian(obji_vec, controls) + SX.zeros(ny, nu);
JxN = jacobian(objN_vec, states) + SX.zeros(nyN, nx);

obji = 0.5*norm_2(obji_vec)^2;
objN = 0.5*norm_2(objN_vec)^2;
gxi = jacobian(obji,states)' + SX.zeros(nx,1);
gui = jacobian(obji,controls)' + SX.zeros(nu,1);
gxN = jacobian(objN,states)' + SX.zeros(nx,1);

Cxi = jacobian(path_con, states) + SX.zeros(nc, nx);
Cui = jacobian(path_con, controls) + SX.zeros(nc, nu);
CxN = jacobian(path_con_N, states) + SX.zeros(ncN, nx);

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

adj_dG = SX.zeros(nx+nu,1) + jtimes(X, [states;controls], lambdai, true);

if nc>0
    adj_dB = SX.zeros(nx+nu,1) + jtimes(path_con, [states;controls], mui, true);
else
    adj_dB = SX.zeros(nx+nu,1);
end

if ncN>0
    adj_dBN = SX.zeros(nx,1) + jtimes(path_con_N, states, muN, true);
else
    adj_dBN = SX.zeros(nx,1);
end

adj_fun = Function('adj_fun',{states,controls,params,refs,Q,lambdai,mui},{dobj, adj_dG, adj_dB});
adjN_fun = Function('adjN_fun',{states,params,refN, QN, muN},{dobjN, adj_dBN});

%% Code generation and Compile

generate=input('Would you like to generate the source code?(y/n)','s');

if strcmp(generate,'y')

    display('                           ');
    display('    Generating source code...');
    
    cd model_src
      
    opts = struct( 'main', false, 'mex' , true ) ; 
    Simulate_system.generate('Simulate_system.c',opts);
    h_fun.generate('h_fun.c',opts);
    path_con_fun.generate('path_con_fun.c',opts);
    path_con_N_fun.generate('path_con_N_fun.c',opts);
    Ji_fun.generate('Ji_fun.c',opts);
    JN_fun.generate('JN_fun.c',opts);
   
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
        
        P.generate();
    cd ../model_src

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
       CC_FLAGS = 'GCC="/usr/bin/gcc-4.9"';
    end
    
    OP_FLAGS='-O';
    PRINT_FLAGS='-silent';
    
    mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, 'path_con_fun.c');
    mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, 'path_con_N_fun.c');
    mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, 'h_fun.c');
    mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, 'Simulate_system.c');
    mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, 'Ji_fun.c');
    mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, 'JN_fun.c');
       
    cd ../mex_core
    Compile_Mex;
    cd ../model_src

display('    Compilation completed!');

end

cd ..
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

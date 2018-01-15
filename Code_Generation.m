clear all;clc;
display('-----------------------------------------------------');
display('This framework is developed by Yutao Chen, DEI, UniPD');
display('-----------------------------------------------------');


%% Insert Model here
settings.model='Hexacopter';

switch settings.model
    case 'InvertedPendulum'
        InvertedPendulum;
    case 'DiM'
        DiM;
    case 'ChainofMasses_Lin'
        ChainofMasses_Lin;
    case 'Hexacopter'
        Hexacopter;
end

%%
import casadi.*

lambdai=SX.sym('lambdai',nx,1);            % the i th multiplier for equality constraints
mui=SX.sym('mui',nc,1);                  % the i th multiplier for inequality constraints
muN=SX.sym('muN',ncN,1);                 % the N th multiplier for inequality constraints

Ts    = 0.01;  % NMPC sampling time [s]

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

Ts_st = 0.1; % shooting interval time
s  = 2; % No. of integration steps per shooting interval
DT = Ts_st/s;
f_fun  = Function('f_fun', {states,controls,params}, {SX.zeros(nx,1)+x_dot},{'states','controls','params'},{'xdot'});
jacX = SX.zeros(nx,nx)+jacobian(x_dot,states);
jacU = SX.zeros(nx,nu)+jacobian(x_dot,controls);
jac_f_fun  = Function('jac_f_fun', {states,controls,params}, {jacX,jacU});

impl_jac_x = SX.zeros(nx,nx)+jacobian(impl_f,states);
impl_jac_u = SX.zeros(nx,nu)+jacobian(impl_f,controls);
impl_jac_xdot = SX.zeros(nx,nx)+jacobian(impl_f,xdot);
impl_f = SX.zeros(nx,1) + impl_f;
impl_f_fun = Function('impl_f_fun',{states,controls,params,xdot},{impl_f, impl_jac_x, impl_jac_u, impl_jac_xdot});

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
z = [states;controls];
F = Function('F', {z,params}, {X + SX.zeros(nx,1)}, {'z','params'}, {'xf'});
A = jacobian(X,states) + SX.zeros(nx,nx);
B = jacobian(X,controls) + SX.zeros(nx,nu);
D = Function('D', {z,params}, {A, B}, {'z','params'}, {'A','B'});

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

Ji_fun=Function('Ji_fun',{z,params,refs,Q},{Jxi,Jui},{'z','params','refs','Q'},{'Jxi','Jui'});
JN_fun=Function('JN_fun',{states,params,refN,QN},{JxN},{'states','params','refN','QN'},{'JxN'});

gi_fun=Function('gi_fun',{z,params,refs,Q},{gxi, gui},{'z','params','refs','Q'},{'gxi','gui'});
gN_fun=Function('gN_fun',{states,params,refN,QN},{gxN},{'states','params','refN','QN'},{'gN'});

Ci_fun=Function('Ci_fun',{z},{Cxi, Cui},{'z'},{'Cxi','Cui'});
CN_fun=Function('CN_fun',{states},{CxN},{'states'},{'CxN'});

dobj = SX.zeros(nx+nu,1) + jacobian(obji,z)';
dobjN = SX.zeros(nx,1) + jacobian(objN,states)';
adj_dG = SX.zeros(nx+nu,1) + jtimes(X, z, lambdai, true);

if nc>0
    adj_dB = SX.zeros(nx+nu,1) + jtimes(path_con, z, mui, true);
else
    adj_dB =[];
end

if ncN>0
    adj_dBN = SX.zeros(nx,1) + jtimes(path_con_N, states, muN, true);
else
    adj_dBN = [];
end

adj_fun = Function('adj_fun',{z,params,refs,Q, lambdai, mui},{dobj, adj_dG, adj_dB});
adjN_fun = Function('adjN_fun',{states,params,refN, QN, muN},{dobjN, adj_dBN});

%% Code generation and Compile

generate=input('Would you like to generate the source code?(y/n)','s');

if strcmp(generate,'y')

    display('                           ');
    display('Generating source code...');

    if exist([pwd,'/Source_Codes'],'dir')~=7
        mkdir([pwd,'/Source_Codes']);
    end
    
    cd Source_Codes
      
    opts = struct( 'main', false, 'mex' , true ) ; 
    Simulate_system.generate('Simulate_system.c',opts);
    h_fun.generate('h_fun.c',opts);
    path_con_fun.generate('path_con_fun.c',opts);
    path_con_N_fun.generate('path_con_N_fun.c',opts);
   
    opts = struct('main',false,'mex',false,'with_header',true);
    cd ../mex_core
        P = CodeGenerator ('casadi_src.c', opts) ;
        P.add(f_fun);
        P.add(jac_f_fun);
        P.add(vdeFun);
        P.add(impl_f_fun);
        P.add(F);
        P.add(D);
        P.add(h_fun);
        P.add(path_con_fun);
        P.add(path_con_N_fun);
        P.add(gi_fun);
        P.add(gN_fun);
        P.add(Ji_fun);
        P.add(JN_fun);
        P.add(Ci_fun);
        P.add(CN_fun);
        P.add(adj_fun);
        P.add(adjN_fun);
        
        P.generate();
    cd ../Source_Codes

display('                           ');
display('Code generation completed!');

end

compile=input('Would you like to compile the source code?(y/n)','s');
if strcmp(compile,'y')
    
    display('                           ');
    display('Compiling...(turn off optimization flag is recommended)');
    op_flag=input('Would you like to turn OFF the optimization flag?(y/n)','s');

    if isempty(op_flag)
        op_flag = 'y';
    end

    if strcmp(op_flag,'n')
        display('Optimization flag turned on');
        display('                           ');
        
       mex -largeArrayDims path_con_fun.c
       mex -largeArrayDims path_con_N_fun.c
       mex -largeArrayDims h_fun.c
       mex -largeArrayDims Simulate_system.c
       
       cd ../mex_core
       Compile_Mex;
       cd ../Source_Codes
              
    end

    if strcmp(op_flag,'y')
        display('Optimization flag turned off');
        display('                            ');
        
        mex -largeArrayDims -g path_con_fun.c
        mex -largeArrayDims -g path_con_N_fun.c
        mex -largeArrayDims -g h_fun.c
        mex -largeArrayDims -g Simulate_system.c
        
        cd ../mex_core
        Compile_Mex
        cd ../Source_Codes

    end

cd ..
display('                           ');
display('Compilation completed!');

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

save('settings','settings');

clear all;

display('                           ');
display('NMPC solver prepared! Enjoy solving...');
display('                           ');

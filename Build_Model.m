function settings = Build_Model(model, buildMex)
%BUILD_MODEL Generate model C code and optionally build MATMPC MEX files.
%   Build_Model('InvertedPendulum') builds the model and saves its settings.
%   Build_Model(model, false) generates source files only.

if nargin < 1
    model = 'InvertedPendulum';
end
if nargin < 2
    buildMex = true;
end
validateattributes(buildMex, {'logical'}, {'scalar'}, mfilename, 'buildMex');
if isempty(which('casadi.SX'))
    error('MATMPC:MissingCasADi', 'Add your CasADi MATLAB installation to the path before building.');
end
if buildMex && isempty(mex.getCompilerConfigurations('C', 'Selected'))
    error('MATMPC:MissingCompiler', 'Configure a supported C compiler with mex -setup C before building.');
end

projectRoot = fileparts(mfilename('fullpath'));
originalFolder = pwd;
originalPath = path;
stage = tempname;
mkdir(stage);
restore = onCleanup(@() cleanBuild(stage, originalFolder, originalPath));
stagedModel = fullfile(stage, 'model_src');
stagedCore = fullfile(stage, 'mex_core');
mkdir(stagedModel);
mkdir(stagedCore);
addpath(fullfile(projectRoot, 'examples'), '-end');
addpath(fullfile(projectRoot, 'mex_core'), '-begin');

settings.model = char(model);
run(settings.model);

%%
import casadi.*

lambda=SX.sym('lambda',nx,1);            % the i th multiplier for equality constraints
mu_u=SX.sym('mu_u',nu,1);                  % the i th multiplier for bounds on controls
mu_x=SX.sym('mu_x',nbx,1);                  % the i th multiplier for bounds on controls
mu_g=SX.sym('mu_g',nc,1);                  % the i th multiplier for bounds on controls
muN_g=SX.sym('muN_g',ncN,1);                 % the N th multiplier for inequality constraints

%% Generate some functions

f_fun  = Function('f_fun', {states,controls,params,alg}, {SX.zeros(nx,1)+x_dot},{'states','controls','params','alg'},{'xdot'});

impl_jac_f_x = SX.zeros(nx,nx)+jacobian(impl_f,states);
impl_jac_f_u = SX.zeros(nx,nu)+jacobian(impl_f,controls);
impl_jac_f_xdot = SX.zeros(nx,nx)+jacobian(impl_f,xdot);

impl_f_fun = Function('impl_f_fun',{states,controls,params,xdot,alg},{SX.zeros(nx,1) + impl_f});
impl_jac_f_x_fun = Function('impl_jac_f_x_fun',{states,controls,params,xdot,alg},{impl_jac_f_x});
impl_jac_f_u_fun = Function('impl_jac_f_u_fun',{states,controls,params,xdot,alg},{impl_jac_f_u});
impl_jac_f_xdot_fun = Function('impl_jac_f_xdot_fun',{states,controls,params,xdot,alg},{impl_jac_f_xdot});

if nz>0
    g_fun  = Function('g_fun', {states,controls,params,xdot,alg}, {SX.zeros(nz,1)+z_fun},{'states','controls','params','xdot','alg'},{'zfun'});

    impl_jac_f_z = SX.zeros(nx,nz)+jacobian(impl_f,alg);
    impl_jac_g_x = SX.zeros(nz,nx)+jacobian(z_fun,states);
    impl_jac_g_u = SX.zeros(nz,nu)+jacobian(z_fun,controls);
    impl_jac_g_z = SX.zeros(nz,nz)+jacobian(z_fun,alg);

    impl_jac_f_z_fun = Function('impl_jac_f_z_fun',{states,controls,params,xdot,alg},{impl_jac_f_z});
    impl_jac_g_x_fun = Function('impl_jac_g_x_fun',{states,controls,params,xdot,alg},{impl_jac_g_x});
    impl_jac_g_u_fun = Function('impl_jac_g_u_fun',{states,controls,params,xdot,alg},{impl_jac_g_u});
    impl_jac_g_z_fun = Function('impl_jac_g_z_fun',{states,controls,params,xdot,alg},{impl_jac_g_z});
else
    g_fun  = Function('g_fun', {states,controls,params,xdot,alg}, {0},{'states','controls','params','xdot','alg'},{'zfun'});

    impl_jac_f_z_fun = Function('impl_jac_f_z_fun',{states,controls,params,xdot,alg},{0});
    impl_jac_g_x_fun = Function('impl_jac_g_x_fun',{states,controls,params,xdot,alg},{0});
    impl_jac_g_u_fun = Function('impl_jac_g_u_fun',{states,controls,params,xdot,alg},{0});
    impl_jac_g_z_fun = Function('impl_jac_g_z_fun',{states,controls,params,xdot,alg},{0});
end

Sx = SX.sym('Sx',nx,nx);
Su = SX.sym('Su',nx,nu);
if nz==0
    vdeX = SX.zeros(nx,nx);
    vdeX = vdeX + jtimes(x_dot,states,Sx);
    vdeU = SX.zeros(nx,nu) + jacobian(x_dot,controls);
    vdeU = vdeU + jtimes(x_dot,states,Su);
    vdeFun = Function('vdeFun',{states,controls,params,Sx,Su,alg},{vdeX,vdeU});

    adjW = SX.zeros(nx+nu,1) + jtimes(x_dot, [states;controls], lambda, true);
    adj_ERK_fun = Function('adj_ERK_fun',{states,controls,params,lambda,alg},{adjW});
else
    vdeFun = Function('vdeFun',{states,controls,params,Sx,Su,alg},{0,0});
    adj_ERK_fun = Function('adj_ERK_fun',{states,controls,params,lambda,alg},{0});
end

%% objective and constraints

h_fun=Function('h_fun', {states,controls,params}, {h},{'states','controls','params'},{'h'});
hN_fun=Function('hN_fun', {states,params}, {hN},{'states','params'},{'hN'});
path_con_fun=Function('path_con_fun', {states,controls,params}, {general_con},{'states','controls','params'},{'general_con'});
path_con_N_fun=Function('path_con_N_fun', {states,params}, {general_con_N},{'states','params'},{'general_con_N'});

gxi = jacobian(obji,states)' + SX.zeros(nx,1);
gui = jacobian(obji,controls)' + SX.zeros(nu,1);
gxN = jacobian(objN,states)' + SX.zeros(nx,1);

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

obji_fun = Function('obji_fun',{states,controls,params,refs,Q},{obji+SX.zeros(1,1)});
objN_fun = Function('objN_fun',{states,params,refN,QN},{objN+SX.zeros(1,1)});

Hz_fun = Function('Hz_fun',{aux,params,refs,Q},{Hz});
HzN_fun = Function('HzN_fun',{auxN,params,refN,QN},{HzN});
Hw = Hz_fun(h,params,refs,Q) + SX.zeros(ny,ny);
HwN = HzN_fun(hN,params,refN,QN)+ SX.zeros(nyN,nyN);

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

%% Generate model and solver sources

if buildMex
    core = fullfile(projectRoot, 'mex_core');
    copyfile(fullfile(core, '*.c'), stagedCore);
    copyfile(fullfile(core, '*.h'), stagedCore);
end

cd(stagedModel);
opts = struct( 'main', false, 'mex' , true ) ;
h_fun.generate('h_fun.c',opts);
f_fun.generate('f_fun.c',opts);
g_fun.generate('g_fun.c',opts);
path_con_fun.generate('path_con_fun.c',opts);
path_con_N_fun.generate('path_con_N_fun.c',opts);

opts = struct('main',false,'mex',false,'with_header',true);
cd(stagedCore);
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
cd(originalFolder);

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

relativeFiles = {
    fullfile('model_src', 'h_fun.c')
    fullfile('model_src', 'f_fun.c')
    fullfile('model_src', 'g_fun.c')
    fullfile('model_src', 'path_con_fun.c')
    fullfile('model_src', 'path_con_N_fun.c')
    fullfile('mex_core', 'casadi_src.c')
    fullfile('mex_core', 'casadi_src.h')
};

if buildMex
    Compile_Mex(stagedCore, stagedModel);
    extension = ['*.', mexext];
    folders = {'model_src', 'mex_core'};
    for i = 1:numel(folders)
        mexFiles = dir(fullfile(stage, folders{i}, extension));
        for j = 1:numel(mexFiles)
            relativeFiles{end + 1, 1} = fullfile(folders{i}, mexFiles(j).name); %#ok<AGROW>
        end
    end
    mkdir(fullfile(stage, 'data'));
    relativeFiles{end + 1, 1} = fullfile('data', 'settings.mat');
    save(fullfile(stage, relativeFiles{end}), 'settings');
end

publishBuild(stage, projectRoot, relativeFiles);

if buildMex
    fprintf('Built MATMPC MEX files for %s.\n', settings.model);
else
    fprintf('Generated C sources for %s.\n', settings.model);
end
end

function publishBuild(stage, projectRoot, files)
lockFolder = fullfile(projectRoot, '.matmpc-build-lock');
[locked, message] = mkdir(lockFolder);
if ~locked || ~isempty(message)
    error('MATMPC:BuildLock', 'Could not acquire the build lock at %s. %s', lockFolder, message);
end
unlock = onCleanup(@() rmdir(lockFolder));

backup = fullfile(stage, 'previous');
mkdir(backup);
replacing = false(size(files));
for i = 1:numel(files)
    destination = fullfile(projectRoot, files{i});
    [~, mexName, extension] = fileparts(destination);
    if strcmp(extension, ['.', mexext])
        clear(mexName);
    end
    destinationType = exist(destination, 'file');
    if destinationType == 7
        error('MATMPC:BuildDestinationIsDirectory', 'Expected an output file at %s.', destination);
    end
    replacing(i) = destinationType ~= 0;
    if replacing(i)
        copyfile(destination, fullfile(backup, num2str(i)));
    end
end

try
    for i = 1:numel(files)
        copyfile(fullfile(stage, files{i}), fullfile(projectRoot, files{i}), 'f');
    end
catch problem
    for j = i:-1:1
        destination = fullfile(projectRoot, files{j});
        try
            if replacing(j)
                copyfile(fullfile(backup, num2str(j)), destination, 'f');
            elseif exist(destination, 'file') ~= 0
                delete(destination);
            end
        catch rollbackProblem
            problem = addCause(problem, rollbackProblem);
        end
    end
    rethrow(problem);
end
end

function cleanBuild(stage, originalFolder, originalPath)
cd(originalFolder);
path(originalPath);
rmdir(stage, 's');
end

% Initialize either supplied inverted-pendulum Simulink example.

MATMPCController.assertBuildAllowed();
clear mex;

projectRoot = fileparts(mfilename('fullpath'));
settingsFile = fullfile(projectRoot, 'data', 'settings.mat');
if exist(settingsFile, 'file') ~= 2
    error('MATMPC:MissingSimulinkSettings', ...
        'Run Build_Model(''InvertedPendulum'') before initializing the Simulink example.');
end
generated = load(settingsFile, 'settings');
settings = generated.settings;
if ~strcmp(settings.model, 'InvertedPendulum')
    error('MATMPC:InvalidSimulinkModel', ...
        'The supplied Simulink models require settings from Build_Model(''InvertedPendulum'').');
end

addpath(projectRoot, fullfile(projectRoot, 'nmpc'), ...
    fullfile(projectRoot, 'model_src'), fullfile(projectRoot, 'mex_core'));
if ismac
    solverPlatform = 'mac';
elseif isunix
    solverPlatform = 'linux';
elseif ispc
    solverPlatform = 'win64';
else
    error('MATMPC:UnsupportedPlatform', 'No bundled qpOASES build for this platform.');
end
addpath(fullfile(projectRoot, 'solver', solverPlatform, 'qpoases'), '-end');
if exist('qpOASES_sequence', 'file') ~= 3
    error('MATMPC:MissingQpOases', ...
        'Add a qpOASES_sequence MEX file matching %s to the MATLAB path.', mexext);
end

Ts = settings.Ts_st;
Ts_st = settings.Ts_st;
nx = settings.nx;
nu = settings.nu;
nz = settings.nz;
ny = settings.ny;
nyN = settings.nyN;
np = settings.np;
nc = settings.nc;
ncN = settings.ncN;
nbu = settings.nbu;
nbx = settings.nbx;
nbu_idx = settings.nbu_idx;
nbx_idx = settings.nbx_idx;

N = 80;
N2 = N / 5;
r = 10;
settings.N = N;
settings.N2 = N2;
settings.r = r;

opt.hessian = 'Gauss_Newton';
opt.levenberg_marquardt = 0;
opt.integrator = 'ERK4';
opt.condensing = 'default_full';
opt.qpsolver = 'qpoases';
opt.hotstart = 'no';
opt.shifting = 'no';
opt.ref_type = 0;
opt.nonuniform_grid = 0;
opt.RTI = 'yes';

% The charts consume these workspace parameters and the plant uses x0.
[initial, ~] = InitData(settings);
x0 = initial.x0;
u0 = initial.u0;
z0 = initial.z0;
para0 = initial.od(:, 1);
x = initial.x;
u = initial.u;
z = initial.z;
para = initial.od;
W = initial.W;
WN = initial.WN;
lb = initial.lb;
ub = initial.ub;
lbx = initial.lbx;
ubx = initial.ubx;
lbu = initial.lbu;
ubu = initial.ubu;

% The MATLAB Function block parameters need nonempty placeholders when disabled.
if isempty(lb)
    lb = 0;
    ub = 0;
end
if isempty(z0)
    z0 = 0;
    z = zeros(1, N);
end

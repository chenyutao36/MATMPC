function Compile_Mex(coreDir, modelDir)
%COMPILE_MEX Build MATMPC MEX files with MATLAB's configured C compiler.
%   Compile_Mex() builds the core from generated sources in this directory.
%   Compile_Mex(coreDir, modelDir) builds both the model and core MEX files
%   using the portable sources in the sibling core directory.

if exist('MATMPCController', 'class') == 8
    MATMPCController.assertBuildAllowed();
end
if nargin < 1
    coreDir = fileparts(mfilename('fullpath'));
end
if nargin > 1
    modelFiles = {'path_con_fun.c', 'path_con_N_fun.c', 'h_fun.c', 'f_fun.c', 'g_fun.c'};
    for i = 1:numel(modelFiles)
        compile(modelDir, modelFiles(i), {});
    end
end

blas = {'-lmwblas'};
lapack = {'-lmwblas', '-lmwlapack'};
common = {'mpc_common.c'};
sharedCore = fullfile(coreDir, '..', 'core');
sharedInclude = ['-I', sharedCore];
laBackend = fullfile(coreDir, 'matmpc_la_mex.c');
integrationKernels = {sharedInclude, fullfile(sharedCore, 'erk_kernel.c'), ...
    fullfile(sharedCore, 'irk_kernel.c'), fullfile(sharedCore, 'matmpc_lu.c'), laBackend};
recoverKernel = {sharedInclude, fullfile(sharedCore, 'recover_kernel.c')};
condenseKernels = {sharedInclude, fullfile(sharedCore, 'condense_affine.c'), ...
    fullfile(sharedCore, 'condense_matrix.c'), ...
    fullfile(coreDir, 'matmpc_la_mex.c')};
condenseMoveBlockedKernels = {sharedInclude, fullfile(sharedCore, 'condense_affine.c'), ...
    fullfile(sharedCore, 'condense_mb.c'), fullfile(coreDir, 'matmpc_la_mex.c')};
gaussNewtonKernel = fullfile(sharedCore, 'qp_gauss_newton.c');
modelSupport = {'casadi_wrapper.c', 'sim.c', 'erk.c', 'irk_mex.c', 'irk_ode.c', ...
    'irk_dae.c', 'casadi_src.c', 'mpc_common.c'};

compile(coreDir, [{'qp_generation.c'}, modelSupport], lapack, integrationKernels{:}, gaussNewtonKernel);
compile(coreDir, {'Condensing.c'}, blas, condenseKernels{:});
compile(coreDir, [{'sparse2full.c'}, common], {});
compile(coreDir, [{'partial_condensing_default.c', 'partial_condensing_routines.c'}, common], blas);
compile(coreDir, {'Recover.c'}, {}, recoverKernel{:});
compile(coreDir, [{'Line_search.c'}, modelSupport], lapack, integrationKernels{:});
compile(coreDir, [{'solution_info.c'}, modelSupport], lapack, integrationKernels{:});
compile(coreDir, [{'Simulate_System.c'}, modelSupport], lapack, integrationKernels{:});

compile(coreDir, {'Condensing_mb.c'}, blas, condenseMoveBlockedKernels{:});
compile(coreDir, [{'qp_generation_ngrid.c'}, modelSupport], lapack, integrationKernels{:});
compile(coreDir, [{'qp_generation_mb.c'}, modelSupport], lapack, integrationKernels{:});
end

function compile(sourceDir, filenames, libraries, varargin)
sources = fullfile(sourceDir, filenames);
mex('-largeArrayDims', '-O', '-silent', '-outdir', sourceDir, sources{:}, varargin{:}, libraries{:});
end

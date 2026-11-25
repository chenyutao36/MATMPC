function Compile_Mex(coreDir, modelDir)
%COMPILE_MEX Build MATMPC MEX files with MATLAB's configured C compiler.
%   Compile_Mex() builds the core from generated sources in this directory.
%   Compile_Mex(coreDir, modelDir) builds both the model and core MEX files.

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
modelSupport = {'casadi_wrapper.c', 'sim.c', 'erk.c', 'irk_ode.c', ...
    'irk_dae.c', 'casadi_src.c', 'mpc_common.c'};

compile(coreDir, [{'qp_generation.c'}, modelSupport], lapack);
compile(coreDir, [{'Condensing.c'}, common], blas);
compile(coreDir, [{'sparse2full.c'}, common], {});
compile(coreDir, [{'partial_condensing_default.c', 'partial_condensing_routines.c'}, common], blas);
compile(coreDir, [{'Recover.c'}, common], blas);
compile(coreDir, [{'Line_search.c'}, modelSupport], lapack);
compile(coreDir, [{'solution_info.c'}, modelSupport], lapack);
compile(coreDir, [{'Simulate_System.c'}, modelSupport], lapack);

compile(coreDir, [{'Condensing_mb.c'}, common], blas);
compile(coreDir, [{'qp_generation_ngrid.c'}, modelSupport], lapack);
compile(coreDir, [{'qp_generation_mb.c'}, modelSupport], lapack);
end

function compile(sourceDir, filenames, libraries)
sources = fullfile(sourceDir, filenames);
mex('-largeArrayDims', '-O', '-silent', '-outdir', sourceDir, sources{:}, libraries{:});
end

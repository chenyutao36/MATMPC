currentFolder = pwd;
objFolder = '/home/yutaochen/Documents/MATLAB/Packages/MATMPC-GITLAB/MATMPC/mex_core';

if ~strcmp(currentFolder, objFolder)
    cd /home/yutaochen/Documents/MATLAB/Packages/MATMPC-GITLAB/MATMPC/mex_core
end

% detect OS type

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
    disp('Platform not supported')
end

%% Configure Compiler

options = '-largeArrayDims';

if OS_WIN
   CC_FLAGS='';
end
if OS_LINUX 
   CC_FLAGS = 'GCC="/usr/bin/gcc-4.9"';
end

OP_FLAGS='';

LIB1 = '-lmwblas';
LIB2 = '-lmwlapack';
% LIB3 = '-lopenblas';
% LIB4 = '-lblasfeo';

LIB1_PATH = '';
LIB2_PATH = '';
% LIB3_PATH = '-L/home/yutaochen/Downloads/OpenBLAS/build/lib';
% LIB4_PATH = '-L/home/yutaochen/Documents/MATLAB/Packages/blasfeo/lib';

HEAD1_PATH = '';
HEAD2_PATH = '';
% HEAD3_PATH = '-I/home/yutaochen/Downloads/OpenBLAS/build/include';
% HEAD4_PATH='-I/home/yutaochen/Documents/MATLAB/Packages/blasfeo/include';

%% 

mex(options, CC_FLAGS, OP_FLAGS, 'Condensing_mex_new.c', LIB1);

% mex(options, OP_FLAGS, CC_FLAGS, HEAD1_PATH, LIB1_PATH, 'qp_generation_mex.c','casadi_wrapper.c','sim_erk.c','f_fun.c','jac_f_fun.c','vdeFun.c','impl_f_fun.c','F.c','D.c','Ji_fun.c','gi_fun.c','ineq_fun.c','Ci_fun.c','JN_fun.c','gN_fun.c','ineqN_fun.c','CN_fun.c','adj_fun.c','adjN_fun.c',LIB1);

mex(options, OP_FLAGS, CC_FLAGS, HEAD1_PATH, LIB1_PATH, 'qp_generation_mex_erk.c','casadi_wrapper.c','sim_erk.c','f_fun.c','jac_f_fun.c','vdeFun.c','impl_f_fun.c','F.c','D.c','Ji_fun.c','gi_fun.c','ineq_fun.c','Ci_fun.c','JN_fun.c','gN_fun.c','ineqN_fun.c','CN_fun.c','adj_fun.c','adjN_fun.c',LIB1);

% mex(options, OP_FLAGS, CC_FLAGS, HEAD1_PATH, LIB1_PATH, 'qp_generation_mex_irk.c','casadi_wrapper.c','sim_irk.c','f_fun.c','jac_f_fun.c','vdeFun.c','impl_f_fun.c','F.c','D.c','Ji_fun.c','gi_fun.c','ineq_fun.c','Ci_fun.c','JN_fun.c','gN_fun.c','ineqN_fun.c','CN_fun.c','adj_fun.c','adjN_fun.c', LIB1, LIB2);

mex(options, OP_FLAGS, CC_FLAGS, HEAD1_PATH, LIB1_PATH,'Recover_mex.c', LIB1);

mex(options, OP_FLAGS, CC_FLAGS, HEAD1_PATH, LIB1_PATH,'solution_info.c','casadi_wrapper.c','sim_erk.c','f_fun.c','jac_f_fun.c','vdeFun.c','impl_f_fun.c','F.c','D.c','Ji_fun.c','gi_fun.c','ineq_fun.c','Ci_fun.c','JN_fun.c','gN_fun.c','ineqN_fun.c','CN_fun.c','adj_fun.c','adjN_fun.c', LIB1, LIB2);

mex(options, OP_FLAGS, CC_FLAGS, HEAD1_PATH, LIB1_PATH,'Line_search_mex.c', LIB1);

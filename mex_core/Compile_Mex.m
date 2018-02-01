%% set path for your computer

%% detect OS type

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
   CC_FLAGS='CXXFLAGS="$CXXFLAGS -Wall"'; % use MinGW not VS studio
end
if OS_LINUX 
   CC_FLAGS = 'GCC="/usr/bin/gcc-4.9"';
end

OP_FLAGS='-O';
PRINT_FLAGS='-silent';

LIB1 = '-lmwblas';
LIB2 = '-lmwlapack';
% LIB3 = '-lblasfeo';
% LIB4 = '-lhpipm';

% LIB5 = '-lqpOASES_e';

LIB1_PATH = '';
LIB2_PATH = '';
% LIB3_PATH = '-L/home/chen/Documents/Packages/blasfeo/lib';
% LIB4_PATH = '-L/home/chen/Documents/Packages/hpipm/lib';
% LIB5_PATH = '-L/home/yutaochen/Documents/MATLAB/Packages/qpOASES_C/qpOASES/build/lib';

HEAD1_PATH = '';
HEAD2_PATH = '';
% HEAD3_PATH = '-I/home/chen/Documents/Packages/blasfeo/include';
% HEAD4_PATH = '-I/home/chen/Documents/Packages/hpipm/include';
% HEAD5_PATH='-I/home/yutaochen/Documents/MATLAB/Packages/qpOASES_C/qpOASES/include';
%% These functions should be all compiled

% mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, HEAD1_PATH, LIB1_PATH, 'qp_generation.c','casadi_wrapper.c','sim.c','casadi_src.c','mpc_common.c',LIB1, LIB2);

% mex(options, CC_FLAGS, OP_FLAGS, PRINT_FLAGS, 'Condensing.c','mpc_common.c', LIB1);

% mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, HEAD1_PATH, LIB1_PATH,'Recover.c', LIB1);

% mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, HEAD1_PATH, LIB1_PATH,'solution_info.c','casadi_wrapper.c','casadi_src.c','sim.c','mpc_common.c', LIB1, LIB2);

% mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, HEAD1_PATH, LIB1_PATH,'Line_search.c','casadi_wrapper.c','casadi_src.c','sim.c','mpc_common.c', LIB1, LIB2);

%% only for testing, don't touch

mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, HEAD1_PATH, LIB1_PATH, 'qp_generation_cmon.c','casadi_wrapper.c','casadi_src.c','mpc_common.c',LIB1);
% mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, HEAD1_PATH, LIB1_PATH, 'adaptive_eta.c',LIB1);

% mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, 'HPIPM_d_solve_ipm2_hard_ocp_qp.c', LIB3_PATH, LIB3, HEAD3_PATH, LIB4_PATH, LIB4, HEAD4_PATH);
% mex GCC='/usr/bin/gcc-4.9' HPIPM_d_solve_ipm2_hard_ocp_qp.c /home/chen/Documents/Packages/hpipm/lib/libhpipm.a /home/chen/Documents/Packages/blasfeo/lib/libblasfeo.a -I/home/chen/Documents/Packages/hpipm/include -I/home/chen/Documents/Packages/blasfeo/include

% mex(options, OP_FLAGS, CC_FLAGS, HEAD4_PATH, LIB4_PATH, LIB4, 'blasfeo.c');

% mex(options, OP_FLAGS, CC_FLAGS, HEAD5_PATH, LIB5_PATH, LIB5, 'qpoases_c.c');

% mex(options, OP_FLAGS, CC_FLAGS, HEAD5_PATH, LIB5_PATH, LIB5, 'qpoases_hotstart.c');
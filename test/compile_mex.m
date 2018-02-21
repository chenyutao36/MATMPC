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

LIB1_PATH = '';
LIB2_PATH = '';

HEAD1_PATH = '';
HEAD2_PATH = '';

mex(options, OP_FLAGS, CC_FLAGS, PRINT_FLAGS, HEAD1_PATH, LIB1_PATH, 'RTI_step.c','RTI_step_common.c','RTI_step_funcs.c','casadi_wrapper.c','casadi_src.c','mpc_common.c',LIB1, LIB2);
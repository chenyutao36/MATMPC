%% detect OS type

% mexfiles = ['RTI_step.c ', ...
%             ' RTI_step_common.c ',' RTI_step_funcs.c ',' casadi_wrapper.c ',' casadi_src.c ',' mpc_common.c ', ...  
%             ' /home/chen/Documents/Packages/OpenBLAS-0.2.20/build/lib/libopenblas_haswellp-r0.2.20.a '
%            ];
%        
% mexcmd = 'mex';
% 
% mexcmd = [mexcmd, ' -DINT64 CFLAGS="\$CFLAGS -std=c99"'];
% 
% mexcmd = [mexcmd, ' -I.. -I/home/chen/Documents/Packages/OpenBLAS-0.2.20/build/include'];
% 
% mexcmd = [mexcmd, ' ', mexfiles];
% 
% mexcmd = [mexcmd, '-lpthread'];
% 
% eval(mexcmd);

%%

mexfiles = ['RTI_step.c ', ...
            ' RTI_step_common.c ',' RTI_step_funcs.c ',' casadi_wrapper.c ',' casadi_src.c ',' mpc_common.c ', ...  
            ' /home/chen/Documents/Packages/OpenBLAS-0.2.20/build/lib/libopenblas_haswellp-r0.2.20.a ', ...
            ' /home/chen/Documents/Packages/QORE/bin/libqpsolver_dense.a ', ...
            ' /home/chen/Documents/Packages/QORE/bin/libkktpack_dense.a ', ...
            ' /home/chen/Documents/Packages/QORE/bin/libqpcore.a '
           ];
       
mexcmd = 'mex';

mexcmd = [mexcmd, ' -DINT64 CFLAGS="\$CFLAGS -std=c99" GCC="/usr/bin/gcc-4.9"'];

mexcmd = [mexcmd, ' -I.. -I/home/chen/Documents/Packages/OpenBLAS-0.2.20/build/include -I/home/chen/Documents/Packages/QORE/ -I/home/chen/Documents/Packages/QORE/QPSOLVER_DENSE/include'];

mexcmd = [mexcmd, ' ', mexfiles];

mexcmd = [mexcmd, '-lpthread -lrt -lblas -L/home/chen/Documents/Packages/QORE/external/blasfeo/lib -lblasfeo'];

eval(mexcmd);
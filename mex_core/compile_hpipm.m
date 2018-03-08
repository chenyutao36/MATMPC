
%% compile hpipm sparse

% mexfiles = ['hpipm_sparse.c ', ...
%             ' /home/chen/Documents/Packages/OpenBLAS-0.2.20/build/lib/libopenblas_haswellp-r0.2.20.a ', ...
%             ' /home/chen/Documents/Packages/hpipm/lib/libhpipm.a ', ...
%             ' /home/chen/Documents/Packages/blasfeo/lib/libblasfeo.a ',...
%            ];

%% compile hpipm dense

mexfiles = ['hpipm_dense.c ', ...
            ' /home/chen/Documents/Packages/hpipm/lib/libhpipm.a ', ...
            ' /home/chen/Documents/Packages/blasfeo/lib/libblasfeo.a ',...
           ];

%%
       
mexcmd = 'mex';

mexcmd = [mexcmd, ' -DINT64 CFLAGS="\$CFLAGS -std=c99" GCC="/usr/bin/gcc-4.9"'];

mexcmd = [mexcmd, ' -I.. -I/home/chen/Documents/Packages/OpenBLAS-0.2.20/build/include -I/home/chen/Documents/Packages/hpipm/include -I/home/chen/Documents/Packages/blasfeo/include'];

mexcmd = [mexcmd, ' ', mexfiles];

eval(mexcmd);

%% hpipm sparse

mexfiles_sp = ['hpipm_sparse.c ', ...
            ' /opt/hpipm/lib/libhpipm.a ', ...
            ' /opt/blasfeo/lib/libblasfeo.a ',...
           ];

%% hpipm dense

mexfiles_d = ['hpipm_dense.c ', ...
            ' /opt/hpipm/lib/libhpipm.a ', ...
            ' /opt/blasfeo/lib/libblasfeo.a ',...
           ];
       
%% hpipm condense
mexfiles_cond = ['condensing_hpipm.c ', ...
            ' /opt/hpipm/lib/libhpipm.a ', ...
            ' /opt/blasfeo/lib/libblasfeo.a ',...
           ];
       
%% blasfeo condense
mexfiles_bcd = ['condensing_blasfeo.c ', ...
            ' /opt/blasfeo/lib/libblasfeo.a ',...
           ];
       
%% hpipm pcond
mexfiles_pcond = ['hpipm_pcond.c ', ...
            ' /opt/hpipm/lib/libhpipm.a ', ...
            ' /opt/blasfeo/lib/libblasfeo.a ',...
           ];

%%
       
mexcmd = 'mex';

mexcmd = [mexcmd, ' -O -DINT64 CFLAGS="\$CFLAGS -std=c99" GCC="/usr/bin/gcc-4.9"'];

mexcmd = [mexcmd, ' -I.. -I/opt/hpipm/include -I/opt/blasfeo/include'];

%%

mexcmd_sp = [mexcmd, ' ', mexfiles_sp];
mexcmd_d = [mexcmd, ' ', mexfiles_d];
mexcmd_cond = [mexcmd, ' ', mexfiles_cond];
mexcmd_bcd = [mexcmd, ' ', mexfiles_bcd];
mexcmd_pcond = [mexcmd, ' ', mexfiles_pcond];

% eval(mexcmd_sp);
% eval(mexcmd_d);
eval(mexcmd_cond);
% eval(mexcmd_bcd);
% eval(mexcmd_pcond);
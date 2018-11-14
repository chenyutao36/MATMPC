
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
    PREFIX = 'C:\msys64';
else
    disp('Platform not supported')
end

%% hpipm sparse

if OS_LINUX
    mexfiles_sp = ['hpipm_sparse.c ', ...
                ' /opt/hpipm/lib/libhpipm.a ', ...
                ' /opt/blasfeo/lib/libblasfeo.a ',...
                   ];
elseif OS_WIN
    mexfiles_sp = ['hpipm_sparse.c ', ...
                   PREFIX,'\opt\hpipm\lib\libhpipm.a ', ...
                   PREFIX,'\opt\blasfeo\lib\libblasfeo.a ',...
                   ];
end
                   
       
%% hpipm pcond

if OS_LINUX
    mexfiles_pcondsol = ['hpipm_pcond.c ', ...
                ' /opt/hpipm/lib/libhpipm.a ', ...
                ' /opt/blasfeo/lib/libblasfeo.a ',...
                   ];
elseif OS_WIN
    mexfiles_pcondsol = ['hpipm_pcond.c ', ...
                   PREFIX,'\opt\hpipm\lib\libhpipm.a ', ...
                   PREFIX,'\opt\blasfeo\lib\libblasfeo.a ',...
                   ];
end

%% blasfeo condensing

if OS_LINUX
    mexfiles_bcond = ['Condensing_Blasfeo.c ', ...
                ' /opt/blasfeo/lib/libblasfeo.a ',...
                   ];
elseif OS_WIN
    mexfiles_bcond = ['Condensing_Blasfeo.c ', ...
                   PREFIX,'\opt\blasfeo\lib\libblasfeo.a ',...
                   ];
end

%% Build mex command
       
mexcmd = 'mex';

if OS_LINUX
    mexcmd = [mexcmd, ' -O -DINT64 CFLAGS="\$CFLAGS -std=c99" GCC="/usr/bin/gcc"'];
    mexcmd = [mexcmd, ' -I.. -I/opt/hpipm/include -I/opt/blasfeo/include'];
elseif OS_WIN
    mexcmd = [mexcmd, ' -O -DINT64 CFLAGS="$CFLAGS -std=c99" '];
    mexcmd = [mexcmd, ' -I.. -I' PREFIX '\opt\hpipm\include -I' PREFIX '\opt\blasfeo\include'];
end

%%

mexcmd_sp = [mexcmd, ' ', mexfiles_sp];
mexcmd_pcondsol = [mexcmd, ' ', mexfiles_pcondsol];
mexcmd_bcond = [mexcmd, ' ', mexfiles_bcond];

eval(mexcmd_sp);
eval(mexcmd_pcondsol);
eval(mexcmd_bcond);
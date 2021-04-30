
%% detect OS type

OS_MAC = 0;
OS_LINUX = 0;
OS_WIN = 0;

%% Important! 
% Please revise the following codes using the correct 
% path for Blasfeo and Hpipm that have been installed in your PC.

if ismac
    OS_MAC = 1;
    PREFIX = '/Users/xxx'; % xxx is the username of your MAC OS
elseif isunix
    OS_LINUX = 1;
elseif ispc
    OS_WIN = 1;
    PREFIX = 'D:\Tools';
else
    disp('Platform not supported')
end

%% hpipm sparse

if OS_LINUX
    mexfiles_sp = ['hpipm_sparse.c ', ...
                ' /opt/hpipm/lib/libhpipm.a ', ...
                ' /opt/blasfeo/lib/libblasfeo.a ',...
                   ]; % BLASFEO and HPIPM are in default installed at /opt/
elseif OS_WIN
    mexfiles_sp = ['hpipm_sparse.c ', ...
                   PREFIX,'\hpipm\lib\libhpipm.a ', ...
                   PREFIX,'\blasfeo\lib\libblasfeo.a ',...
                   ]; % BLASFEO and HPIPM are in default installed at D:\Tools\
elseif OS_MAC
    mexfiles_sp = ['hpipm_sparse.c ', ...
                    PREFIX,'/Documents/blasfeo/lib/libblasfeo.a ', ...
                    PREFIX,'/Documents/hpipm/lib/libhpipm.a ',...
                       ]; % BLASFEO and HPIPM are in default installed at Users/xxx/Documents/, where xxx is the username of your MAC OS
end
                   
       
%% hpipm pcond

if OS_LINUX
    mexfiles_pcondsol = ['hpipm_pcond.c ', ...
                ' /opt/hpipm/lib/libhpipm.a ', ...
                ' /opt/blasfeo/lib/libblasfeo.a ',...
                   ];
elseif OS_WIN
    mexfiles_pcondsol = ['hpipm_pcond.c ', ...
                   PREFIX,'\hpipm\lib\libhpipm.a ', ...
                   PREFIX,'\blasfeo\lib\libblasfeo.a ',...
                   ];
elseif OS_MAC
    mexfiles_pcondsol = ['hpipm_pcond.c ', ...
                    PREFIX,'/Documents/blasfeo/lib/libblasfeo.a ', ...
                    PREFIX,'/Documents/hpipm/lib/libhpipm.a ',...
                       ];
end

%% blasfeo condensing
% 
% if OS_LINUX
%     mexfiles_bcond = ['Condensing_Blasfeo.c ', ...
%                       ' /opt/blasfeo/lib/libblasfeo.a ',...
%                       ];
% elseif OS_WIN
%     mexfiles_bcond = ['Condensing_Blasfeo.c ', ...
%                       PREFIX,'\opt\blasfeo\lib\libblasfeo.a ',...
%                       ];
% elseif OS_MAC
%     mexfiles_bcond = ['Condensing_Blasfeo.c ', ...
%                       PREFIX,'/Documents/blasfeo_lib/blasfeo/lib/libblasfeo.a ',...
%                       ];
% end

%% Build mex command
       
mexcmd = 'mex';

if OS_LINUX
    mexcmd = [mexcmd, ' -O -DINT64 CFLAGS="\$CFLAGS -std=c99" GCC="/usr/bin/gcc"'];
    mexcmd = [mexcmd, ' -I.. -I/opt/hpipm/include -I/opt/blasfeo/include'];
elseif OS_WIN
    mexcmd = [mexcmd, ' -O -DINT64 CFLAGS="$CFLAGS -std=c99" '];
    mexcmd = [mexcmd, ' -I.. -I' PREFIX '\hpipm\include -I' PREFIX '\blasfeo\include'];
elseif OS_MAC
    mexcmd = [mexcmd, ' -O -DINT64 CFLAGS="\$CFLAGS -std=c99"'];
    mexcmd = [mexcmd, ' -I.. -I' PREFIX, '/Documents/hpipm/include -I' PREFIX '/Documents/blasfeo/include'];
end

%%

mexcmd_sp = [mexcmd, ' ', mexfiles_sp];
mexcmd_pcondsol = [mexcmd, ' ', mexfiles_pcondsol];
% mexcmd_bcond = [mexcmd, ' ', mexfiles_bcond];

eval(mexcmd_sp);
eval(mexcmd_pcondsol);
% eval(mexcmd_bcond);
%% Initialization
function [input, data] = InitData(settings)

    Ts  = settings.Ts;       % Sampling time
    Ts_st = settings.Ts_st;  % Shooting interval
    s = settings.s;      % number of integration steps per interval
    nx = settings.nx;    % No. of states
    nu = settings.nu;    % No. of controls
    ny = settings.ny;    % No. of outputs (references)    
    nyN= settings.nyN;   % No. of outputs at terminal stage 
    np = settings.np;    % No. of parameters (on-line data)
    nc = settings.nc;    % No. of constraints
    ncN = settings.ncN;  % No. of constraints at terminal stage
    N     = settings.N;             % No. of shooting points
    nbx = settings.nbx;
    nbu = settings.nbu;
    nbu_idx = settings.nbu_idx;

    switch settings.model
        case 'DiM'
            input.x0 = zeros(nx,1);    % initial state
            input.u0 = zeros(nu,1);    % initial control
            para0 = 0;  % initial parameters (by default a np by 1 vector, if there is no parameter, set para0=0)

            %weighting matrices
            Q=diag([1200,1200,2000,800,800,5800,... % perceived acc and angular vel
                    32000*1.1,32000*1.1,1600*1,... %px,py,pz hex
                    3200*1.1,3200*1.1,2000*1,... %vx, vy, vz hex
                    4600*1,600*1,... % x,y tri
                    850*1,850*1,... % vx,vy tri
                    3700,3000,1500,... % phi, theta, psi hex
                    750,... % phi tri
                    0.01,0.0,0.0,... % omega phi,theta,psi hex
                    500.0,... % omega phi tri
                    0.0,0.0,0.001,... %ax,ay,az hex %         20*1.1,20*1.1,... % ax,ay tri
                    0.0,0.01,0.1 ... % alpha phi,theta, psi hex 
                    ]);

              QN=Q(1:nyN,1:nyN);
              
              % upper and lower bounds for states (=nbx)
              lb_x = [];
              ub_x = [];
              lb_xN = [];
              ub_xN = [];
              
              % upper and lower bounds for controls (=nbu)           
              lb_u = [];
              ub_u = [];

              % upper and lower bounds for general constraints (=nc)
              lb_g=[1.045;1.045;1.045;1.045;1.045;1.045];    % lower bounds for ineq constraints
              ub_g=[1.3750;1.3750;1.3750;1.3750;1.3750;1.3750];  % upper bounds for ineq constraints
              lb_gN=[1.045;1.045;1.045;1.045;1.045;1.045];  % lower bounds for ineq constraints at terminal point
              ub_gN=[1.3750;1.3750;1.3750;1.3750;1.3750;1.3750];  % upper bounds for ineq constraints at terminal point
              
              % store the constraint data into input
              input.lb=repmat([lb_g;lb_x],1,N);
              input.ub=repmat([ub_g;ub_x],1,N); 
              input.lbN=[lb_gN;lb_xN];               
              input.ubN=[ub_gN;ub_xN]; 

              lbu = -inf(nu,1);
              ubu = inf(nu,1);
              for i=1:nbu
                  lbu(nbu_idx(i)) = lb_u(i);
                  ubu(nbu_idx(i)) = ub_u(i);
              end

              input.lbu = repmat(lbu,1,N);
              input.ubu = repmat(ubu,1,N);

        case 'InvertedPendulum'
            input.x0 = [0;pi;0;0];    
            input.u0 = zeros(nu,1);    
            para0 = 0;  

            Q=diag([10 10 0.1 0.1 0.01]);
            QN=Q(1:nyN,1:nyN);

            % upper and lower bounds for states (=nbx)
            lb_x = -2;
            ub_x = 2;
            lb_xN = -2;
            ub_xN = 2;

            % upper and lower bounds for controls (=nbu)           
            lb_u = -20;
            ub_u = 20;
                       
            % upper and lower bounds for general constraints (=nc)
            lb_g = [];
            ub_g = [];            
            lb_gN = [];
            ub_gN = [];

            % store the constraint data into input
            input.lb=repmat([lb_g;lb_x],1,N);
            input.ub=repmat([ub_g;ub_x],1,N); 
            input.lbN=[lb_gN;lb_xN];               
            input.ubN=[ub_gN;ub_xN]; 
            
            lbu = -inf(nu,1);
            ubu = inf(nu,1);
            for i=1:nbu
                lbu(nbu_idx(i)) = lb_u(i);
                ubu(nbu_idx(i)) = ub_u(i);
            end
            
            input.lbu = repmat(lbu,1,N);
            input.ubu = repmat(ubu,1,N);

        case 'ChainofMasses_Lin'
            n=5;
            data.n=n;
            input.x0=zeros(nx,1);
            for i=1:n
                input.x0(i)=7.5*i/n;
            end
            input.u0=zeros(nu,1);
            para0=0;
            wv=[];wx=[];wu=[];
            wu=blkdiag(wu,0.1, 0.1, 0.1);
            for i=1:3
                wx=blkdiag(wx,25);
                wv=blkdiag(wv,diag(0.25*ones(1,n-1)));
            end
            Q=blkdiag(wx,wv,wu);
            QN=blkdiag(wx,wv);

            % upper and lower bounds for states (=nbx)
            lb_x = [];
            ub_x = [];
            lb_xN = [];
            ub_xN = [];

            % upper and lower bounds for controls (=nbu)           
            lb_u = [-1;-1;-1];
            ub_u = [1;1;1];
                       
            % upper and lower bounds for general constraints (=nc)
            lb_g = [];
            ub_g = [];            
            lb_gN = [];
            ub_gN = [];

            % store the constraint data into input
            input.lb=repmat([lb_g;lb_x],1,N);
            input.ub=repmat([ub_g;ub_x],1,N); 
            input.lbN=[lb_gN;lb_xN];               
            input.ubN=[ub_gN;ub_xN]; 
            
            lbu = -inf(nu,1);
            ubu = inf(nu,1);
            for i=1:nbu
                lbu(nbu_idx(i)) = lb_u(i);
                ubu(nbu_idx(i)) = ub_u(i);
            end
            
            input.lbu = repmat(lbu,1,N);
            input.ubu = repmat(ubu,1,N);

        case 'ChainofMasses_NLin'
            n=10;
            data.n=n;
    %         x0=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 zeros(1,nx-n)]';
            input.x0=[rand(1,n), 0.6*rand(1,n)-1, -0.6*rand(1,n) , zeros(1,3*(n-1))]';
    %         xref=[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 zeros(1,nx-n)]';
            input.u0=zeros(nu,1);
            para0=0;
            wv=[];wx=[];wu=[];
            wu=blkdiag(wu,0.01, 0.01, 0.01);
            for i=1:3
                wx=blkdiag(wx,25);
                wv=blkdiag(wv,diag(1*ones(1,n-1)));
            end
            Q=blkdiag(wx,wv,wu);
            QN=blkdiag(wx,wv);

            % upper and lower bounds for states (=nbx)
            lb_x = [];
            ub_x = [];
            lb_xN = [];
            ub_xN = [];

            % upper and lower bounds for controls (=nbu)           
            lb_u = [-1;-1;-1];
            ub_u = [1;1;1];
                       
            % upper and lower bounds for general constraints (=nc)
            lb_g = [];
            ub_g = [];            
            lb_gN = [];
            ub_gN = [];

            % store the constraint data into input
            input.lb=repmat([lb_g;lb_x],1,N);
            input.ub=repmat([ub_g;ub_x],1,N); 
            input.lbN=[lb_gN;lb_xN];               
            input.ubN=[ub_gN;ub_xN]; 
            
            lbu = -inf(nu,1);
            ubu = inf(nu,1);
            for i=1:nbu
                lbu(nbu_idx(i)) = lb_u(i);
                ubu(nbu_idx(i)) = ub_u(i);
            end
            
            input.lbu = repmat(lbu,1,N);
            input.ubu = repmat(ubu,1,N);

        case 'Hexacopter'
            input.x0=zeros(nx,1);
            input.u0=zeros(nu,1);
            para0=0;
            
            q = [5e0, 5e0, 5e0, 0.1, 0.1, 0.1];
            qN = q(1:nyN);
            Q = diag(q);
            QN = diag(qN);

            % upper and lower bounds for states (=nbx)
            lb_x = [];
            ub_x = [];

            % upper and lower bounds for controls (=nbu)           
            lb_u = -inf*ones(nbu,1);
            ub_u = inf*ones(nbu,1);
                       
            % upper and lower bounds for general constraints (=nc)
            lb_g = [];
            ub_g = [];            
            lb_gN = [];
            ub_gN = [];

            % store the constraint data into input
            input.lb=repmat([lb_g;lb_x],1,N);
            input.ub=repmat([ub_g;ub_x],1,N); 
            input.lbN=[lb_gN;lb_x];               
            input.ubN=[ub_gN;ub_x]; 
            
            lbu = -inf(nu,1);
            ubu = inf(nu,1);
            for i=1:nbu
                lbu(nbu_idx(i)) = lb_u(i);
                ubu(nbu_idx(i)) = ub_u(i);
            end
            
            input.lbu = repmat(lbu,1,N);
            input.ubu = repmat(ubu,1,N);

        case 'TiltHex'
            input.x0=zeros(nx,1);
            input.u0=zeros(nu,1);
            para0=0;

            q =[5,5,5,0.1,1,0.1,1e-5*ones(1,nu)];
            qN = q(1:nyN);
            Q = diag(q);
            QN = diag(qN);
            
            % upper and lower bounds for states (=nbx)
            lb_x = 0*ones(nbx,1);
            ub_x = 12*ones(nbx,1);

            % upper and lower bounds for controls (=nbu)           
            lb_u = -80*ones(nbu,1);
            ub_u = 80*ones(nbu,1);
                       
            % upper and lower bounds for general constraints (=nc)
            lb_g = [];
            ub_g = [];            
            lb_gN = [];
            ub_gN = [];

            % store the constraint data into input
            input.lb=repmat([lb_g;lb_x],1,N);
            input.ub=repmat([ub_g;ub_x],1,N); 
            input.lbN=[lb_gN;lb_x];               
            input.ubN=[ub_gN;ub_x]; 
            
            lbu = -inf(nu,1);
            ubu = inf(nu,1);
            for i=1:nbu
                lbu(nbu_idx(i)) = lb_u(i);
                ubu(nbu_idx(i)) = ub_u(i);
            end
            
            input.lbu = repmat(lbu,1,N);
            input.ubu = repmat(ubu,1,N);

            %Frequency for x(t) in rad/s 
            data.f_rif_x=1.2;
            %Same frequency used in MPC algorithm
            data.f_x=data.f_rif_x*0.5/pi;
            %Amplitude of x(t)
            data.amplitude_x=1.2;
            %Frequency for theta(t) in rad/s
            data.f_rif_theta=1.2;
            %Same frequency used in MPC algorithm
            data.f_theta=data.f_rif_theta*0.5/pi;
            %Amplitude of theta(t)
            data.amplitude_theta=pi/18;
    end

    % prepare the data

    x = repmat(input.x0,1,N+1);  % initialize all shooting points with the same initial state 
    % x = repmat(xref,1,N+1); 
    u = repmat(input.u0,1,N);    % initialize all controls with the same initial control
    para = repmat(para0,1,N+1); % initialize all parameters with the same initial para

    input.z=[x(:,1:N);u];        % states and controls of the first N stages (N by (nx+nu) matrix)
    input.xN=x(:,N+1);          % states of the terminal stage (nx by 1 vector)
    input.od=para;               % on-line parameters (N+1 by np matrix)
    input.W=Q;                   % weights of the first N stages (ny by ny matrix)
    input.WN=QN;                 % weights of the terminal stage (nyN by nyN matrix)


    %% Reference generation

    switch settings.model
        case 'DiM'

            load REF_DiM_2;

            REF_DiM_2 = [REF_DiM_2, zeros(5000,24)];

            data.REF = REF_DiM_2;

        case 'InvertedPendulum'

            data.REF=zeros(1,nx+nu);
            
%             T = 5/Ts;
%             data.REF = [zeros(T,1), pi*ones(T,1), zeros(T,3);
%                         1.5*ones(T,1), zeros(T,4);
%                         -1.5*ones(T,1), pi*ones(T,1), zeros(T,3);
%                         1.5*ones(T,1),zeros(T,4);
%                         zeros(10*T,1), pi*ones(10*T,1), zeros(10*T,3)];

        case 'ChainofMasses_Lin'

            data.REF=[7.5,0,0,zeros(1,3*(n-1)),zeros(1,nu)];

        case 'ChainofMasses_NLin'

            data.REF=[1,0,0,zeros(1,3*(n-1)),zeros(1,nu)];

        case 'Hexacopter'

            data.REF = [1 1 1 0 0 0];

        case 'TiltHex'
            data.REF = [];
        
    end
    
end
   
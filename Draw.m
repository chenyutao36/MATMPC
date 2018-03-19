%% Plot results

% define colors
blue=[0 0 255]/255;
red=[220 20 60]/255;
orange=[255 165 0]/255;
green=[0 205 102]/255;

% start generating pictures
switch settings.model
    case 'DiM'

        samples=size(y_sim,1);

        figure;
        title('MCA Tracking of perceived signals');

        subplot(3,2,1);
        plot(y_sim(:,1),'r');
        hold on;
        plot(data.REF(1:samples,1),'k');
        title('Longitudinal: $\hat{a}_x$','Interpreter','latex');

        subplot(3,2,2);
        plot(y_sim(:,2),'r');
        hold on;
        plot(data.REF(1:samples,2),'k');
        title('Lateral: $\hat{a}_y$','Interpreter','latex');

        subplot(3,2,3);
        plot(y_sim(:,3),'r');
        hold on;
        plot(data.REF(1:samples,3),'k');
        title('Vertical: $\hat{a}_z$','Interpreter','latex');

        subplot(3,2,4);
        plot(y_sim(:,4),'r');
        hold on;
        plot(data.REF(1:samples,4),'k');
        title('Roll: $\hat{\omega}_{\psi}$','Interpreter','latex');

        subplot(3,2,5);
        plot(y_sim(:,5),'r');
        hold on;
        plot(data.REF(1:samples,5),'k');
        title('Pitch: $\hat{\omega}_{\theta}$','Interpreter','latex');

        subplot(3,2,6);
        plot(y_sim(:,6),'r');
        hold on;
        plot(data.REF(1:samples,6),'k');
        title('Yaw: $\hat{\omega}_{\phi}$','Interpreter','latex');


        figure;
        subplot(3,2,1)
        plot(y_sim(:,13));
        hold on;
        plot(y_sim(:,7),'r');
        plot(zeros(Tf*100,1),'k');
        title('Longitudinal displacement')
        lgd=legend('tripod: $p_{x,T}$','hex: $p_{x,H}$','ref: $p_{x,T}$');
        set(lgd,'Interpreter','latex');

        subplot(3,2,2)
        plot(y_sim(:,14));
        hold on;
        plot(y_sim(:,8),'r');
        plot(zeros(Tf*100,1),'k');
        title('Lateral displacement')
        lgd=legend('tripod: $p_{y,T}$','hex: $p_{y,H}$','ref: $p_{y,T}$');
        set(lgd,'Interpreter','latex');

        subplot(3,2,3)
        plot(y_sim(:,9),'r');
        hold on;
        plot(zeros(Tf*100,1),'k');
        title('Vertical displacement: $p_{z,H}$','Interpreter','latex');

        subplot(3,2,4)
        plot(y_sim(:,20));
        hold on;
        plot(y_sim(:,17),'r');
        plot(zeros(Tf*100,1),'k');
        title('Yaw')
        lgd=legend('tripod: $\phi_T$','hex: $\phi_H$','ref');
        set(lgd,'Interpreter','latex');

        subplot(3,2,5)
        plot(y_sim(:,18),'r');
        hold on;
        plot(zeros(Tf*100,1),'k');
        title('Pitch');
        lgd=legend('hexpod: $\theta_H$','ref');
        set(lgd,'Interpreter','latex');

        subplot(3,2,6)
        plot(y_sim(:,19),'r');
        hold on;
        plot(zeros(Tf*100,1),'k');
        title('Roll');
        lgd=legend('hexpod: $\psi_H$','ref');
        set(lgd,'Interpreter','latex');

        figure;
        title('Hex actuator constraints')
        plot(constraints(:,1:6));
        hold on;
        plot(1.045*ones(samples,1),':');
        plot(1.375*ones(samples,1),':');
        axis([0 mem.iter 1.0 1.4]);
        title('Hexpod actuator constraints');

        % figure;
        % plot(constraints(:,7:9));
        % hold on;
        % plot(2.5*ones(samples,1),':');
        % plot(4*ones(samples,1),':');
        % title('Tripod actuator constraints');

    case 'InvertedPendulum'
        
        figure(1);
        subplot(321)
        plot(time,state_sim(:,1));
        title('p');
        subplot(322)
        plot(time,state_sim(:,2)*180/pi);
        title('\theta');
        subplot(323)
        plot(time,state_sim(:,3));
        title('v');
        subplot(324)
        plot(time,state_sim(:,4)*180/pi);
        title('\omega');
        subplot(3,2,[5 6]);
        title('F');
        stairs(time,controls_MPC(:,1));
        
    case 'ChainofMasses_Lin'
        figure(1);
        subplot(311);
        plot(time,controls_MPC(:,1));
        ylabel('$u_x$','Interpreter','latex');
        subplot(312);
        plot(time,controls_MPC(:,2));
        ylabel('$u_y$','Interpreter','latex');
        subplot(313);
        plot(time,controls_MPC(:,3));
        ylabel('$u_z$','Interpreter','latex');

        n = data.n;
        figure(2);
        plot3([0,state_sim(1,1:n)], [0,state_sim(1,n+1:2*n)], [0,state_sim(1,2*n+1:3*n)],'Color',red,'LineStyle','--');
        hold on;
        grid on;
        plot3([0,state_sim(end,1:n)], [0,state_sim(end,n+1:2*n)],[0,state_sim(end,2*n+1:3*n)],'Color',blue,'LineStyle','-');
        scatter3([0,state_sim(1,1:n)], [0,state_sim(1,n+1:2*n)], [0,state_sim(1,2*n+1:3*n)],10,'MarkerFaceColor','none');
        scatter3([0,state_sim(end,1:n)], [0,state_sim(end,n+1:2*n)],[0,state_sim(end,2*n+1:3*n)],10,'MarkerFaceColor',red);
        xlabel('X[m]');
        ylabel('Y[m]');
        zlabel('Z[m]');
        
        xlim([-0.5 8]);
        ylim([-1 1]);
        zlim([-6 0.5]);
        
   case 'ChainofMasses_NLin'
        figure(1);
        subplot(311);
        plot(time,controls_MPC(:,1));
        ylabel('$u_x$','Interpreter','latex');
        subplot(312);
        plot(time,controls_MPC(:,2));
        ylabel('$u_y$','Interpreter','latex');
        subplot(313);
        plot(time,controls_MPC(:,3));
        ylabel('$u_z$','Interpreter','latex');

        n = data.n;
        figure(2);
        plot3([0,state_sim(1,1:n)], [0,state_sim(1,n+1:2*n)], [0,state_sim(1,2*n+1:3*n)],'Color',red,'LineStyle','--');
        hold on;
        grid on;
        plot3([0,state_sim(end,1:n)], [0,state_sim(end,n+1:2*n)],[0,state_sim(end,2*n+1:3*n)],'Color',blue,'LineStyle','-');
        scatter3([0,state_sim(1,1:n)], [0,state_sim(1,n+1:2*n)], [0,state_sim(1,2*n+1:3*n)],10,'MarkerFaceColor','none');
        scatter3([0,state_sim(end,1:n)], [0,state_sim(end,n+1:2*n)],[0,state_sim(end,2*n+1:3*n)],10,'MarkerFaceColor',red);
        xlabel('X[m]');
        ylabel('Y[m]');
        zlabel('Z[m]');
        
        xlim([-1.2 1.2]);
        ylim([-1.2 1.2]);
%         zlim([-1.2 1.2]);
        
    case 'Hexacopter'
        figure();      
        subplot(321)
        hold on;
        grid on;
        plot(time(1:end-1),y_sim(:,1),'Color',red);
        title('x');
        
        subplot(322)
        hold on;
        grid on;
        plot(time(1:end-1),y_sim(:,2),'Color',red);
        if isempty(ref_traj) ~=1
            plot(time(1:end-1),ref_traj, 'k--');
        end
        title('y');
        
        subplot(323)
        hold on;
        grid on;
        plot(time(1:end-1),y_sim(:,3),'Color',red);
        title('z');
        
        subplot(324)
        hold on;
        grid on;
        plot(time(1:end-1),y_sim(:,4),'Color',red);
        title('phi');
        
        subplot(325)
        hold on;
        grid on;
        plot(time(1:end-1),y_sim(:,5),'Color',red);
        title('theta');
        
        subplot(326)
        hold on;
        grid on;
        plot(time(1:end-1),y_sim(:,6),'Color',red);
        title('psi');
        
        figure();
        subplot(321)
        hold on;
        grid on;
        plot(time(1:end),controls_MPC(:,1),'Color',red);
        title('f1');
        
        subplot(322)
        hold on;
        grid on;
        plot(time(1:end),controls_MPC(:,2),'Color',red);
        title('f2');
        
        subplot(323)
        hold on;
        grid on;
        plot(time(1:end),controls_MPC(:,3),'Color',red);
        title('f3');
        
        subplot(324)
        hold on;
        grid on;
        plot(time(1:end),controls_MPC(:,4),'Color',red);
        title('f4');
        
        subplot(325)
        hold on;
        grid on;
        plot(time(1:end),controls_MPC(:,5),'Color',red);
        title('f5');
        
        subplot(326)
        hold on;
        grid on;
        plot(time(1:end),controls_MPC(:,6),'Color',red);
        title('f6');
        
        
%         figure();
%         subplot(321)
%         hold on;
%         grid on;
%         plot(time(1:end),state_sim(:,13),'Color',red);
%         title('f1');
%         
%         subplot(322)
%         hold on;
%         grid on;
%         plot(time(1:end),state_sim(:,14),'Color',red);
%         title('f2');
%         
%         subplot(323)
%         hold on;
%         grid on;
%         plot(time(1:end),state_sim(:,15),'Color',red);
%         title('f3');
%         
%         subplot(324)
%         hold on;
%         grid on;
%         plot(time(1:end),state_sim(:,16),'Color',red);
%         title('f4');
%         
%         subplot(325)
%         hold on;
%         grid on;
%         plot(time(1:end),state_sim(:,17),'Color',red);
%         title('f5');
%         
%         subplot(326)
%         hold on;
%         grid on;
%         plot(time(1:end),state_sim(:,18),'Color',red);
%         title('f6');
        
%         figure();
%         subplot(321)
%         hold on;
%         grid on;
%         plot(time(1:end-1),y_sim(:,7),'Color',red);
%         title('df1');
%         
%         subplot(322)
%         hold on;
%         grid on;
%         plot(time(1:end-1),y_sim(:,8),'Color',red);
%         title('df2');
%         
%         subplot(323)
%         hold on;
%         grid on;
%         plot(time(1:end-1),y_sim(:,9),'Color',red);
%         title('df3');
%         
%         subplot(324)
%         hold on;
%         grid on;
%         plot(time(1:end-1),y_sim(:,10),'Color',red);
%         title('df4');
%         
%         subplot(325)
%         hold on;
%         grid on;
%         plot(time(1:end-1),y_sim(:,11),'Color',red);
%         title('df5');
%         
%         subplot(326)
%         hold on;
%         grid on;
%         plot(time(1:end-1),y_sim(:,12),'Color',red);
%         title('df6');

    case 'TiltHex'
        figure();      
        subplot(321)
        hold on;
        grid on;
        plot(time(1:end-1),y_sim(:,1),'Color',red);
        if isempty(ref_traj) ~=1
            plot(time(1:end-1),ref_traj(1,:), 'k--');
        end
        title('x');
        legend('x','ref');
        
        subplot(322)
        hold on;
        grid on;
        plot(time(1:end-1),y_sim(:,2),'Color',red);     
        title('y');
        
        subplot(323)
        hold on;
        grid on;
        plot(time(1:end-1),y_sim(:,3),'Color',red);
        title('z');
        
        subplot(324)
        hold on;
        grid on;
        plot(time(1:end-1),y_sim(:,4),'Color',red);
        title('phi');
        
        subplot(325)
        hold on;
        grid on;
        plot(time(1:end-1),y_sim(:,5),'Color',red);
        if isempty(ref_traj) ~=1
            plot(time(1:end-1),ref_traj(5,:), 'k--');
        end
        title('theta');
        legend('theta','ref');
        
        subplot(326)
        hold on;
        grid on;
        plot(time(1:end-1),y_sim(:,6),'Color',red);
        title('psi');
        
        figure();
        subplot(321)
        hold on;
        grid on;
        plot(time(1:end),state_sim(:,13),'Color',red);
        title('f1');
        
        subplot(322)
        hold on;
        grid on;
        plot(time(1:end),state_sim(:,14),'Color',red);
        title('f2');
        
        subplot(323)
        hold on;
        grid on;
        plot(time(1:end),state_sim(:,15),'Color',red);
        title('f3');
        
        subplot(324)
        hold on;
        grid on;
        plot(time(1:end),state_sim(:,16),'Color',red);
        title('f4');
        
        subplot(325)
        hold on;
        grid on;
        plot(time(1:end),state_sim(:,17),'Color',red);
        title('f5');
        
        subplot(326)
        hold on;
        grid on;
        plot(time(1:end),state_sim(:,18),'Color',red);
        title('f6');
        
        figure();
        subplot(321)
        hold on;
        grid on;
        plot(time(1:end),controls_MPC(:,1),'Color',red);
        title('df1');
        
        subplot(322)
        hold on;
        grid on;
        plot(time(1:end),controls_MPC(:,2),'Color',red);
        title('df2');
        
        subplot(323)
        hold on;
        grid on;
        plot(time(1:end),controls_MPC(:,3),'Color',red);
        title('df3');
        
        subplot(324)
        hold on;
        grid on;
        plot(time(1:end),controls_MPC(:,4),'Color',red);
        title('df4');
        
        subplot(325)
        hold on;
        grid on;
        plot(time(1:end),controls_MPC(:,5),'Color',red);
        title('df5');
        
        subplot(326)
        hold on;
        grid on;
        plot(time(1:end),controls_MPC(:,6),'Color',red);
        title('df6');
        
    case 'ActiveSeat'
        % load the data you saved
        %load('e:/study/NNID/MATMPC/data/ActiveSeat/AS_REF_DATA');
        
        load('C:\Users\enrico\Documents\MATLAB\GITLAB\MATMPC\data\ActiveSeat\AS_REF_DATA');
        figure
        plot(time(1:end-1),y_sim(:,1))
        hold on; grid on;
        plot(time(1:end-1),REF_XY(1,:))
        legend('Actual pitch perceived vel.','Reference pitch perceived vel.');
        xlabel('time [s]'); ylabel('[rad/s]');
        ylim([-0.2 0.2])

        figure
        plot(time(1:end-1),y_sim(:,2))
        hold on; grid on;
        plot(time(1:end-1),REF_XY(2,:))
        legend('Actual longitudinal perceived acc.','Reference longitudinal perceived acc.')
        xlabel('time [s]'); ylabel('[m/s^2]');

        figure
        plot(time(1:end-1),y_sim(:,3))
        hold on; grid on;
        plot(time(1:end-1),REF_XY(3,:))
        legend('Actual pitch displacement','Reference pitch displacement')
        xlabel('time [s]'); ylabel('[rad]');

        figure
        plot(time(1:end-1),y_sim(:,4))
        hold on; grid on;
        plot(time(1:end-1),REF_XY(4,:))
        plot(time(1:end-1),0.75*ones(1,Tf/Ts),'k--')
        plot(time(1:end-1),-0.75*ones(1,Tf/Ts),'k--')
        legend('Actual longitudinal displacement','Reference longitudinal displacement')
        xlabel('time [s]'); ylabel('[m]');
        ylim([-0.8 0.8])

        figure
        plot(time(1:end-1),y_sim(:,5))
        grid on;
        legend('Actual longitudinal velocity')
        xlabel('time [s]'); ylabel('[m/s]');
        ylim([-0.5 0.5])

        % plot y-roll

        figure
        plot(time(1:end-1),y_sim(:,7))
        hold on; grid on;
        plot(time(1:end-1),REF_YX(1,:))
        legend('Actual roll perceived vel.','Reference roll perceived vel.');
        xlabel('time [s]'); ylabel('[rad/s]');
        ylim([-0.2 0.2])

        figure
        plot(time(1:end-1),y_sim(:,8))
        hold on; grid on;
        plot(time(1:end-1),REF_YX(2,:))
        legend('Actual lateral perceived acc.','Reference lateral perceived acc.');
        xlabel('time [s]'); ylabel('[m/s^2]');
        ylim([-1 1])

        figure
        plot(time(1:end-1),y_sim(:,9))
        hold on; grid on;
        plot(time(1:end-1),REF_YX(3,:))
        legend('Actual roll displacement','Reference roll displacement')
        xlabel('time [s]'); ylabel('[rad]');
        ylim([-0.02 0.02])

        figure
        plot(time(1:end-1),y_sim(:,10))
        hold on; grid on;
        plot(time(1:end-1),REF_YX(4,:))
        plot(time(1:end-1),0.75*ones(1,Tf/Ts),'k--')
        plot(time(1:end-1),-0.75*ones(1,Tf/Ts),'k--')
        legend('Actual lateral displacement','Reference lateral displacement')
        xlabel('time [s]'); ylabel('[m]');
        ylim([-0.8 0.8])

        figure
        plot(time(1:end-1),y_sim(:,11))
        grid on;
        legend('Actual lateral velocity')
        xlabel('time [s]'); ylabel('[m/s]');
        
        figure;
        plot(time(1:end-1),rif_pressione(1:Tf/0.005)/0.016) % 0.016 � l'area del cuscinetto
        hold on; grid on;
        plot(time(1:end),input_u(:,7)/0.016);
        plot(time(1:end-1),y_sim(:,19)/0.016);
        legend('Actual lateral trunk pressure','Active seat trunk pressure','Pressure induced by platform motion')
        xlabel('time [s]'); ylabel('[Pa]');
        
%     case 'ActiveSeat_onlyAS'
%         % load the data you saved
%         load('C:\Users\enrico\Documents\MATLAB\GITLAB\MATMPC\data\ActiveSeat_vp\AS_REF_DATA');
%         
%         figure(1)
%         plot(time(1:end-1),y_sim(:,1)/0.016,'r','Linewidth',2)
%         hold on; grid on;
%         plot(time(1:end-1),rif_pressione(1:Tf/0.005)/0.016,'b-.','Linewidth',2)
%         legend('Active seat trunk pressure','Reference Pressure')
%         xlabel('time [s]'); ylabel('[Pa]');
%         xlabel('time [s]'); ylabel('[m/s^2]');
% 
%         figure(2)
%         plot(time(1:end-1),y_sim(:,2),'r','Linewidth',2)
%         hold on; grid on;
%         plot(time(1:end-1),accX(:),'b-.','Linewidth',2)
%         legend('Actual acceleration X','Reference acceleration X')
%         xlabel('time [s]'); ylabel('[m/s^2]');
% 
%         figure(3)
%         plot(time(1:end-1),y_sim(:,3),'r','Linewidth',2)
%         hold on; grid on;
%         plot(time(1:end-1),roll_ref(:),'b-.','Linewidth',2)
%         legend('Actual roll angle','Reference roll angle')
%         xlabel('time [s]'); ylabel('[rad]');
%         
%         figure(4)
%         plot(time(1:end-1),y_sim(:,4),'r','Linewidth',2)
%         grid on;hold on;
%         plot(time(1:end-1),accY(:),'b-.','Linewidth',2)
%         legend('Actual acceleration Y','Reference acceleration Y')
%         xlabel('time [s]'); ylabel('[m/s^2]');
%         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   case 'ActiveSeat_onlyP'
        % load the data you saved
        load([pwd,'/data/ActiveSeat_onlyP/AS_REF_DATA_onlyP']);
        
        figure(1)
        plot(time(1:end-1),y_sim(:,1)/0.016,'r','Linewidth',2)
        hold on; grid on;
        plot(time(1:end-1),rif_pressione(1:Tf/0.005)/0.016,'b-.','Linewidth',2)
        legend('Active seat trunk pressure','Reference Pressure')
        xlabel('time [s]'); ylabel('[Pa]');
        xlabel('time [s]'); ylabel('[m/s^2]');

    case 'TethUAV'
        axes_ref = [];
        axes_lim = [];
        
        figure();      
        subplot(221)
        hold on;
        grid on;
        plot(time(1:end-1),y_sim(:,1),'Color',red);
        if isempty(ref_traj) ~=1
            plot(time(1:end-1),ref_traj(1,:), 'k--');
        end
        title('\phi');
        legend('\phi','ref');
        ax = gca; % current axes
        axes_ref = [axes_ref; ax];
        axes_lim = [axes_lim; ax.YLim];
        
        subplot(222)
        hold on;
        grid on;
        plot(time(1:end-1),y_sim(:,2),'Color',red);
        if isempty(ref_traj) ~=1
            plot(time(1:end-1),ref_traj(2,:), 'k--');
        end
        title('\phi_{dot}');
        legend('\phi_{dot}','ref');
        ax = gca; % current axes
        axes_ref = [axes_ref; ax];
        axes_lim = [axes_lim; ax.YLim];
         
        subplot(223)
        hold on;
        grid on;
        plot(time(1:end-1),y_sim(:,3),'Color',red);
        if isempty(ref_traj) ~=1
            plot(time(1:end-1),ref_traj(3,:), 'k--');
        end
        title('\theta');
        legend('\theta','ref');
        ax = gca; % current axes
        axes_ref = [axes_ref; ax];
        axes_lim = [axes_lim; ax.YLim];
        
        subplot(224)
        hold on;
        grid on;
        plot(time(1:end-1),y_sim(:,4),'Color',red);
        if isempty(ref_traj) ~=1
            plot(time(1:end-1),ref_traj(4,:), 'k--');
        end
        title('\theta_{dot}');
        legend('\theta_{dot}','ref');
        ax = gca; % current axes
        axes_ref = [axes_ref; ax];
        axes_lim = [axes_lim; ax.YLim];
        
        % set axes limits, the same for all the plots:
        maxY = max(axes_lim(:,2));
        minY = min(axes_lim(:,1));
        for i = 1 : length(axes_lim)
            cur_ax = axes_ref(i);
            cur_ax.YLim = [minY maxY];
        end
        
        figure()
        subplot(211)
        hold on;
        grid on;
        plot(time(1:end),state_sim(:,5),'Color',red);
        title('f1');
        
        subplot(212)
        hold on;
        grid on;
        plot(time(1:end),state_sim(:,6),'Color',red);
        title('f2');
        
        figure();
        subplot(211)
        hold on;
        grid on;
        plot(time(1:end),controls_MPC(:,1),'Color',red);
        title('df1');
        
        subplot(212)
        hold on;
        grid on;
        plot(time(1:end),controls_MPC(:,2),'Color',red);
        title('df2');
        
        figure();
        grid on;
        plot(time(1:end-1),constraints(:,1),'Color',red);
        title('fL');
        
    case 'TethUAV_param'
         
        phi_ref = input.od(1,1);
        phi_ref = repmat(phi_ref, size(time));
        theta_ref = input.od(2,1);
        theta_ref = repmat(theta_ref, size(time));
        axes_ref = [];
        axes_lim = [];
        
        figure();      
        subplot(221)
        hold on;
        grid on;
        plot(time(1:end),rad2deg(state_sim(:,1)),'Color',red);
%         if isempty(ref_traj) ~=1
%             plot(time(1:end-1),ref_traj(1,:), 'k--');
%         end
        plot(time(1:end),rad2deg(phi_ref),'k--');
        title('\phi');
        legend('\phi','ref');
        ax = gca; % current axes
        axes_ref = [axes_ref; ax];
        axes_lim = [axes_lim; ax.YLim];
        
        subplot(222)
        hold on;
        grid on;
        plot(time(1:end),rad2deg(state_sim(:,2)),'Color',red);
%         if isempty(ref_traj) ~=1
%             plot(time(1:end-1),ref_traj(2,:), 'k--');
%         end
        title('\phi_{dot}');
        legend('\phi_{dot}')
%         legend('\phi_{dot}','ref');
        ax = gca; % current axes
        axes_ref = [axes_ref; ax];
        axes_lim = [axes_lim; ax.YLim];
         
        subplot(223)
        hold on;
        grid on;
        plot(time(1:end),rad2deg(state_sim(:,3)),'Color',red);
%         if isempty(ref_traj) ~=1
%             plot(time(1:end-1),ref_traj(3,:), 'k--');
%         end
        plot(time(1:end),rad2deg(theta_ref),'k--');
        title('\theta');
        legend('\theta','ref');
        ax = gca; % current axes
        axes_ref = [axes_ref; ax];
        axes_lim = [axes_lim; ax.YLim];
        
        subplot(224)
        hold on;
        grid on;
        plot(time(1:end),rad2deg(state_sim(:,4)),'Color',red);
%         if isempty(ref_traj) ~=1
%             plot(time(1:end-1),ref_traj(4,:), 'k--');
%         end
        title('\theta_{dot}');
        legend('\theta_{dot}');%,'ref');
        ax = gca; % current axes
        axes_ref = [axes_ref; ax];
        axes_lim = [axes_lim; ax.YLim];
        
        % set axes limits, the same for all the plots:
        maxY = max(axes_lim(:,2));
        minY = min(axes_lim(:,1));
        for i = 1 : length(axes_lim)
            cur_ax = axes_ref(i);
            cur_ax.YLim = [minY maxY];
        end
        
        figure()
        subplot(211)
        hold on;
        grid on;
        plot(time(1:end),state_sim(:,5),'Color',red);
        title('f1');
        
        subplot(212)
        hold on;
        grid on;
        plot(time(1:end),state_sim(:,6),'Color',red);
        title('f2');
        
        figure();
        subplot(211)
        hold on;
        grid on;
        plot(time(1:end),controls_MPC(:,1),'Color',red);
        title('df1');
        
        subplot(212)
        hold on;
        grid on;
        plot(time(1:end),controls_MPC(:,2),'Color',red);
        title('df2');
        
        figure();
        grid on;
        plot(time(1:end-1),constraints(:,1),'Color',red);
        title('fL');

        % plot time statistics
        
        figure();
        hold on;
        grid on;
        plot(time(2:end-1)', CPT(2:end, 1)); % cpt, tshooting, tcond, tqp
%         plot(time(2:end-1)', CPT(2:end, 2));
%         plot(time(2:end-1)', CPT(2:end, 3));
%         plot(time(2:end-1)', CPT(2:end, 4));
        title('Time statistics');
%         legend('cpt', 'tshooting', 'tcond', 'tqp');
        xlabel('[s]')
        ylabel('[ms]')
        
        % plot nonlinear cost fcn terms
        
        figure();
        hold on;
        grid on;
        plot(time', rad2deg(state_sim(:,1)+state_sim(:,3)));
        plot(time', rad2deg(pi/2)*ones(size(time)));
        title('Cost fcn: Avoid singularity');
        legend('\phi + \theta', '\pi/2');
        xlabel('[s]')
        ylabel('[deg]')
        
        figure();
        hold on;
        grid on;
        plot(time', rad2deg(state_sim(:,1)), 'k');
        plot(time, rad2deg(phi_ref), 'k--');
        plot(time', rad2deg(state_sim(:,3)), 'r');
        plot(time, rad2deg(theta_ref), 'r--');
        title('Cost fcn: Attitude behavior close to the ground');
        legend('\phi', '\phi_{ref}', '\theta', '\theta_{ref}');
        xlabel('[s]')
        ylabel('[deg]')
        
        figure();
        hold on;
        grid on;
        plot(time', rad2deg(state_sim(:,2)), 'r');
        plot(time', rad2deg(state_sim(:,4)), 'b');
        plot(time', rad2deg(state_sim(:,1)), 'k');
        plot(time, rad2deg(phi_ref), 'k--');
        title('Cost fcn: \phi and \theta velocities close to the ground');
        leg = legend('$\dot{\phi}$', '$\dot{\theta}$', '$\phi$', '$\phi_{ref}$');
        set(leg,'Interpreter','latex');
        xlabel('[s]')
        ylabel('[deg]')
        
    case 'TethUAV_param_1order_slack'
         
        phi_ref = input.od(1,1);
        phi_ref = repmat(phi_ref, size(time));
        theta_ref = input.od(2,1);
        theta_ref = repmat(theta_ref, size(time));
        axes_ref = [];
        axes_lim = [];
        
        figure();      
        subplot(221)
        hold on;
        grid on;
        plot(time(1:end),rad2deg(state_sim(:,1)),'Color',red);
%         if isempty(ref_traj) ~=1
%             plot(time(1:end-1),ref_traj(1,:), 'k--');
%         end
        plot(time(1:end),rad2deg(phi_ref),'k--');
        title('\phi');
        legend('\phi','ref');
        ax = gca; % current axes
        axes_ref = [axes_ref; ax];
        axes_lim = [axes_lim; ax.YLim];
        
        subplot(222)
        hold on;
        grid on;
        plot(time(1:end),rad2deg(state_sim(:,2)),'Color',red);
%         if isempty(ref_traj) ~=1
%             plot(time(1:end-1),ref_traj(2,:), 'k--');
%         end
        title('\phi_{dot}');
        legend('\phi_{dot}')
%         legend('\phi_{dot}','ref');
        ax = gca; % current axes
        axes_ref = [axes_ref; ax];
        axes_lim = [axes_lim; ax.YLim];
         
        subplot(223)
        hold on;
        grid on;
        plot(time(1:end),rad2deg(state_sim(:,3)),'Color',red);
%         if isempty(ref_traj) ~=1
%             plot(time(1:end-1),ref_traj(3,:), 'k--');
%         end
        plot(time(1:end),rad2deg(theta_ref),'k--');
        title('\theta');
        legend('\theta','ref');
        ax = gca; % current axes
        axes_ref = [axes_ref; ax];
        axes_lim = [axes_lim; ax.YLim];
        
        subplot(224)
        hold on;
        grid on;
        plot(time(1:end),rad2deg(state_sim(:,4)),'Color',red);
%         if isempty(ref_traj) ~=1
%             plot(time(1:end-1),ref_traj(4,:), 'k--');
%         end
        title('\theta_{dot}');
        legend('\theta_{dot}');%,'ref');
        ax = gca; % current axes
        axes_ref = [axes_ref; ax];
        axes_lim = [axes_lim; ax.YLim];
        
        % set axes limits, the same for all the plots:
        maxY = max(axes_lim(:,2));
        minY = min(axes_lim(:,1));
        for i = 1 : length(axes_lim)
            cur_ax = axes_ref(i);
            cur_ax.YLim = [minY maxY];
        end
        
        figure()
        subplot(211)
        hold on;
        grid on;
        plot(time(1:end),state_sim(:,5),'Color',red);
        title('f1');
        
        subplot(212)
        hold on;
        grid on;
        plot(time(1:end),state_sim(:,6),'Color',red);
        title('f2');
        
        figure();
        subplot(211)
        hold on;
        grid on;
        plot(time(1:end),controls_MPC(:,1),'Color',red);
        title('df1');
        
        subplot(212)
        hold on;
        grid on;
        plot(time(1:end),controls_MPC(:,2),'Color',red);
        title('df2');
        
        figure();
        grid on;
        plot(time(1:end-1),constraints(:,1),'Color',red);
        title('fL');

        % plot time statistics
        
        figure();
        hold on;
        grid on;
        plot(time(2:end-1)', CPT(2:end, 1)); % cpt, tshooting, tcond, tqp
%         plot(time(2:end-1)', CPT(2:end, 2));
%         plot(time(2:end-1)', CPT(2:end, 3));
%         plot(time(2:end-1)', CPT(2:end, 4));
        title('Time statistics');
%         legend('cpt', 'tshooting', 'tcond', 'tqp');
        xlabel('[s]')
        ylabel('[ms]')
        
        % plot nonlinear cost fcn terms
        
        figure();
        hold on;
        grid on;
        plot(time', rad2deg(state_sim(:,1)+state_sim(:,3)));
        plot(time', rad2deg(pi/2)*ones(size(time)));
        title('Cost fcn: Avoid singularity');
        legend('\phi + \theta', '\pi/2');
        xlabel('[s]')
        ylabel('[deg]')
        
        figure();
        hold on;
        grid on;
        plot(time', rad2deg(state_sim(:,1)), 'k');
        plot(time, rad2deg(phi_ref), 'k--');
        plot(time', rad2deg(state_sim(:,3)), 'r');
        plot(time, rad2deg(theta_ref), 'r--');
        title('Cost fcn: Attitude behavior close to the ground');
        legend('\phi', '\phi_{ref}', '\theta', '\theta_{ref}');
        xlabel('[s]')
        ylabel('[deg]')
        
        figure();
        hold on;
        grid on;
        plot(time', rad2deg(state_sim(:,2)), 'r');
        plot(time', rad2deg(state_sim(:,4)), 'b');
        plot(time', rad2deg(state_sim(:,1)), 'k');
        plot(time, rad2deg(phi_ref), 'k--');
        title('Cost fcn: \phi and \theta velocities close to the ground');
        leg = legend('$\dot{\phi}$', '$\dot{\theta}$', '$\phi$', '$\phi_{ref}$');
        set(leg,'Interpreter','latex');
        xlabel('[s]')
        ylabel('[deg]')

end
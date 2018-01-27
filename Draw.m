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
        plot(REF(1:samples,1),'k');
        title('Longitudinal: $\hat{a}_x$','Interpreter','latex');

        subplot(3,2,2);
        plot(y_sim(:,2),'r');
        hold on;
        plot(REF(1:samples,2),'k');
        title('Lateral: $\hat{a}_y$','Interpreter','latex');

        subplot(3,2,3);
        plot(y_sim(:,3),'r');
        hold on;
        plot(REF(1:samples,3),'k');
        title('Vertical: $\hat{a}_z$','Interpreter','latex');

        subplot(3,2,4);
        plot(y_sim(:,4),'r');
        hold on;
        plot(REF(1:samples,4),'k');
        title('Roll: $\hat{\omega}_{\psi}$','Interpreter','latex');

        subplot(3,2,5);
        plot(y_sim(:,5),'r');
        hold on;
        plot(REF(1:samples,5),'k');
        title('Pitch: $\hat{\omega}_{\theta}$','Interpreter','latex');

        subplot(3,2,6);
        plot(y_sim(:,6),'r');
        hold on;
        plot(REF(1:samples,6),'k');
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
        axis([0 iter 1.0 1.4]);
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
end


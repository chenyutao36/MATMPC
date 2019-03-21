%% Plot results

% define colors
blue=[0 0 255]/255;
red=[220 20 60]/255;
orange=[255 165 0]/255;
green=[0 205 102]/255;

% start generating pictures
switch settings.model
    
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
                
    case 'TethUAV'
         
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
        title('\phi_{dot}');
        legend('\phi_{dot}')
        ax = gca; % current axes
        axes_ref = [axes_ref; ax];
        axes_lim = [axes_lim; ax.YLim];
         
        subplot(223)
        hold on;
        grid on;
        plot(time(1:end),rad2deg(state_sim(:,3)),'Color',red);
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
        title('Time statistics');
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
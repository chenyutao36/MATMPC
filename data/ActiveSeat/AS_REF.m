function [Yref] = AS_REF(Tf,Ts)

    %% find your path to the original active seat model files
    cd('e:/study/NNID/Active_seat_belt/Nonlinear');

    %% read data
    filen = 'Calabogie';
    sdata = resreader('C_Segment_Calabogie_Drive_Close_200Hz.res');
    
    GenerateInitialValueForCoder_test;
    
    global Var_xy_ Var_yx_ Var_za_ Var_zv_ status_xy status_yx status_za status_zv;
    Var_xy_ = Var_xy;
    Var_yx_ = Var_yx;
    Var_za_ = Var_za;
    Var_zv_ = Var_zv;
    status_xy = 0; 
    status_yx = 0;
    status_za = 0;
    status_zv = 0;


%% Input signals assignment
  
    time = sdata.time_TIME;
    chassis_accelerations_longitudinal = sdata.chassis_accelerations_longitudinal; 
    chassis_accelerations_lateral = sdata.chassis_accelerations_lateral;
    chassis_accelerations_vertical = sdata.chassis_accelerations_vertical;
    chassis_velocities_yaw = sdata.chassis_velocities_yaw;
    chassis_displacements_pitch = sdata.chassis_displacements_pitch;
    chassis_displacements_roll = sdata.chassis_displacements_roll;
    chassis_displacements_roll_wrt_road = sdata.chassis_displacements_roll_wrt_road;
    chassis_displacements_pitch_wrt_road = sdata.chassis_displacements_pitch_wrt_road;
    chassis_displacements_vertical_wrt_road = sdata.chassis_displacements_vertical_wrt_road;

    IN1_XY = interp1(time,chassis_accelerations_longitudinal,time(1):Var_xy.Ts:time(end),'spline');
    IN1_YX = -interp1(time,chassis_accelerations_lateral,time(1):Var_xy.Ts:time(end),'spline');
    IN1_ZV = interp1(time,chassis_accelerations_vertical,time(1):Var_xy.Ts:time(end),'spline');

    IN1_ZA = interp1(time,chassis_velocities_yaw,time(1):Var_xy.Ts:time(end),'spline');
    IN2_XY = interp1(time,chassis_displacements_pitch(1:end),time(1):Var_xy.Ts:time(end),'spline');
    IN2_YX = -interp1(time,chassis_displacements_roll(1:end),time(1):Var_xy.Ts:time(end),'spline');
    
    %% generate reference
    REF_XY = zeros(4,Tf*100); REFP_XY = zeros(2,Tf*100); y_xy = zeros(6,Tf*100); u_xy = zeros(2,Tf*100);
    REF_YX = zeros(4,Tf*100); REFP_YX = zeros(2,Tf*100); y_yx = zeros(6,Tf*100); u_yx = zeros(2,Tf*100);
    REF_ZA = zeros(2,Tf*100); REFP_ZA = zeros(1,Tf*100); y_za = zeros(3,Tf*100); u_za = zeros(1,Tf*100);
    REF_ZV = zeros(2,Tf*100); REFP_ZV = zeros(1,Tf*100); y_zv = zeros(3,Tf*100); u_zv = zeros(1,Tf*100);
    
    load rif_pressione_calabogie.mat
    
    Yref=[];
    for i=1:Tf/Ts
        in_xy = [IN1_XY(i),IN2_XY(i),i];
        [rif_corr_xy, rif_tot_x, rif_perc_x] = funzMC_XY(in_xy); 
        
        REF_XY(:,i) = rif_tot_x(1:4,:); % gli ultimi 2 sono 0 per def
        REFP_XY(:,i) = rif_perc_x(1:2,:); % mi interessano solo vel.ang e acc
        
        in_yx = [IN1_YX(i),IN2_YX(i),i];
        [rif_corr_yx, rif_tot_y, rif_perc_y] = funzMC_YX(in_yx);    

        REF_YX(:,i) = rif_tot_y(1:4,:); % gli ultimi 2 sono 0 per def
        REFP_YX(:,i) = rif_perc_y(1:2,:); % mi interessano solo vel.ang e acc
        
        in_za = [IN1_ZA(i),i];
        [rif_perc_za, rif_tot_za] = funzMC_ZA(in_za);   

        REF_ZA(:,i) = rif_tot_za(1:2,:); % l'ultimo è 0 per def
        REFP_ZA(:,i) = rif_perc_za(1,:); % mi interessa solo vel.ang
        
        in_zv = [IN1_ZV(i),i];
        [rif_perc_zv, rif_tot_zv] = funzMC_ZV(in_zv);  

        REF_ZV(:,i) = rif_tot_zv(1:2,:); % l'ultimo è 0 per def
        REFP_ZV(:,i) = rif_perc_zv(1,:); % mi interessa solo acc
        
        Yref = [Yref; rif_corr_xy' rif_corr_yx' rif_tot_zv' rif_tot_za' rif_pressione(i), zeros(1,7), zeros(1,6)];
    end

    %% save your data in the path of your MATMPC
    save('e:/study/NNID/MATMPC/data/ActiveSeat/AS_REF_DATA', 'REF_XY', 'REF_YX', 'rif_pressione');
    
    cd('e:/study/NNID/MATMPC');
    
    clc;
end


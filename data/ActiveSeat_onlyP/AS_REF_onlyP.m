function [Yref] = AS_REF_onlyP(Tf,Ts)

    %% find your path to the original active seat model files
%     cd(['C:\Users\enrico\Documents\MATLAB\GITLAB\MATMPC\data\ActiveSeat_onlyP']);

% %% Input signals assignment
    
    load rif_pressione_calabogie.mat
    Yref=[];


for i=1:Tf/Ts
            Yref = [Yref; rif_pressione(i),0];
            
end

    %% save your data in the path of your MATMPC
%     save(['C:\Users\enrico\Documents\MATLAB\GITLAB\MATMPC\data\ActiveSeat_onlyP\AS_REF_DATA_onlyP'], 'rif_pressione');
    
%     cd('C:\Users\enrico\Documents\MATLAB\GITLAB\MATMPC'); %return to main
    
%     clc;
end


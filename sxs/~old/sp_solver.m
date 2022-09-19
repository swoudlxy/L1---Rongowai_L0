
%function [sx_pos_xyz,inc_angle_deg,dis_to_coast_km] = sp_solver(tx,rx,ddm,dem1,dem2,dtu10,dist_to_coast_nz)
clear
clc

% load DTU10 model
dtu_path = '../dat/dtu/';
dtu_filename = 'dtu10_v1.dat';

dtu10 = get_dtu10([dtu_path dtu_filename]);

load('../exp/tx1.mat');
load('../exp/rx1.mat');
load('../exp/ddm1.mat');

tx = tx1;   rx = rx1;   ddm = ddm1;

clear tx1 rx1 ddm1

%%
clc

% retrieve tx and rx position vectors
tx_pos_xyz = tx.tx_pos_xyz;
rx_pos_xyz = rx.rx_pos_xyz;

% step 0 - check if LOS exists
LOS_flag = los(tx_pos_xyz,rx_pos_xyz);      % LOS flag, 1 = LOS exit

if LOS_flag == 1
        
    % step 1 - derive inital SP coordinate on WGS84
    [~,sx_lla_coarse] = coarsetune(tx_pos_xyz,rx_pos_xyz);    

    % step 2 - fine-tune initial SP coordinate
    L_ocean_deg = 1;                        % search area starting from 1*1 deg
    res_ocean_meter = 0.01;                 % targeted grid resolution - 1 cm

    [sx_pos_xyz,inc_angle_deg] = finetune_ocean(tx_pos_xyz,rx_pos_xyz, ...
        sx_lla_coarse,dtu10,L_ocean_deg,res_ocean_meter);

    sx_pos_lla = ecef2lla(sx_pos_xyz);

elseif LOS_flag == 0

    sx_pos_xyz = [-9999,-9999,-9999];
    inc_angle_deg = -9999;
    dis_to_coast_km = -9999;

end

%% validation
clc


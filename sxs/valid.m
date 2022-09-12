clear
clc

load('../int/tx_pos_xyz')
load('../int/tx_vel_xyz')

load('../int/prn_code')

load('../int/rx_pos_xyz')
load('../int/rx_vel_xyz')

load('../int/sp_pos_xyz')

load('../int/gps_off_boresight_angle_deg')

load('../int/sp_theta_orbit')
load('../int/sp_az_orbit')

load('../int/sp_theta_body')
load('../int/sp_az_body')

load('../int/static_gps_eirp')

gps_path = '../dat/gps/';
SV_PRN_filename = 'PRN_SV_LUT_v1.csv';
SV_eirp_filename = 'GPS_SV_EIRP_Params_v7.csv';

SV_eirp_LUT = readmatrix([gps_path SV_eirp_filename]);
SV_PRN_LUT = readmatrix([gps_path SV_PRN_filename]);
SV_PRN_LUT = SV_PRN_LUT(:,1:2);

I = length(gps_off_boresight_angle_deg);

for i = 1:4000%I

    tx_pos_xyz1 = tx_pos_xyz(i,:);
    tx_vel_xyz1 = tx_vel_xyz(i,:);

    PRN = prn_code(i);
    
    rx_pos_xyz1 = rx_pos_xyz(i,:);
    rx_vel_xyz1 = rx_vel_xyz(i,:);

    sx_pos_xyz1 = sp_pos_xyz(i,:);

    [gps_off_boresight_angle1,~] = ecef2orf(tx_pos_xyz1,tx_vel_xyz1,sx_pos_xyz1);
    [theta1_orf,phi1_orf] = ecef2orf(rx_pos_xyz1,rx_vel_xyz1,sx_pos_xyz1);
    [stat_eirp_watt1,~,~] = eirp(PRN,gps_off_boresight_angle1,SV_PRN_LUT,SV_eirp_LUT);

    gps_off_boresight_angle(i) = gps_off_boresight_angle1;
    theta_orf(i) = theta1_orf;
    phi_orf(i) = phi1_orf;
    stat_eirp_watt(i) = stat_eirp_watt1;

end

d1 = gps_off_boresight_angle_deg(1:4000)-gps_off_boresight_angle';
d2 = sp_theta_orbit(1:4000)-theta_orf';
d3 = sp_az_orbit(1:4000)-phi_orf';
d4 = static_gps_eirp(1:4000)-stat_eirp_watt';
%%
close all

figure
subplot(2,2,1),plot(d1),grid on
xlabel('Samples'),ylim([-0.02 0.02])
title('Difference - Off-boresight Angle in Degrees')

subplot(2,2,2),plot(d4),grid on
xlabel('Samples'),ylim([-0.02 0.02])
title('Difference - Static GPS EIRP in Watts')

subplot(2,2,3),plot(d2),grid on
xlabel('Samples'),ylim([-0.02 0.02])
title('Difference - Ele Angle along SP in SC BOF in Degrees')

subplot(2,2,4),plot(d3),grid on
xlabel('Samples'),ylim([-0.02 0.02])
title('Difference - Az Angle along SP in SC BOF in Degrees')
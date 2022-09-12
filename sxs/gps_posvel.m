% this function generates the position and velocity vectors in ECEF
% coordinate frame for a given SV's PRN and a given GPS timestamp
% this function also generates 

function [tx_pos_xyz,tx_vel_xyz,gps_clk_bias,gps_clk_drift] = gps_posvel(PRN,gps_timestamp,sp3_filename)

gps_week = gps_timestamp.gps_week;
gps_tow = gps_timestamp.gps_tow;

tx_info = GPS_GetSVInfo(PRN,gps_week,gps_tow,sp3_filename);

tx_pos_xyz = tx_info(1:3);
tx_vel_xyz = tx_info(5:7);

gps_clk_bias = tx_info(4);
gps_clk_drift = tx_info(8);
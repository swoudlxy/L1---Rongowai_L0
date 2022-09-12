% This function interpolates values outputed by NGRx 
% stamped by PVT epoch to ddm epoch
% Inputs:
% 1) pvt_utc,ddm_utc: pvt and ddm epoch
% 2) rx_pos_xyz_pvt,rx_vel_xyz_pvt,rx_attitude_pvt,rx_clk_pvt:
% rx position, velocity, attitude, and clock vectors, 
% stamped at PVT epoch
% Outputs:
% 1) rx_pos_xyz,rx_vel_xyz,rx_attitude,rx_clk:
% rx position, velocity, attitude, and clock vectors,
% stamped at ddm epoch

function [rx_pos_xyz,rx_vel_xyz,rx_attitude,rx_clk] = PVT2ddm_timestamp(pvt_utc,ddm_utc,rx_pos_xyz_pvt, ...
    rx_vel_xyz_pvt,rx_attitude_pvt,rx_clk_pvt)

rx_pos_x_pvt = rx_pos_xyz_pvt(:,1);
rx_pos_y_pvt = rx_pos_xyz_pvt(:,2);
rx_pos_z_pvt = rx_pos_xyz_pvt(:,3);

rx_vel_x_pvt = rx_vel_xyz_pvt(:,1);
rx_vel_y_pvt = rx_vel_xyz_pvt(:,2);
rx_vel_z_pvt = rx_vel_xyz_pvt(:,3);

rx_roll_deg_pvt = rx_attitude_pvt(:,1);
rx_pitch_deg_pvt = rx_attitude_pvt(:,2);
rx_yaw_deg_pvt = rx_attitude_pvt(:,3);

rx_clk_bias_m_pvt = rx_clk_pvt(:,1);
rx_clk_drift_mps_pvt = rx_clk_pvt(:,2);

% linear interpolation all the values at ddm timestamp
rx_pos_x = interp1(pvt_utc,rx_pos_x_pvt,ddm_utc);
rx_pos_y = interp1(pvt_utc,rx_pos_y_pvt,ddm_utc);
rx_pos_z = interp1(pvt_utc,rx_pos_z_pvt,ddm_utc);

rx_pos_xyz = [rx_pos_x,rx_pos_y,rx_pos_z];

rx_vel_x = interp1(pvt_utc,rx_vel_x_pvt,ddm_utc);
rx_vel_y = interp1(pvt_utc,rx_vel_y_pvt,ddm_utc);
rx_vel_z = interp1(pvt_utc,rx_vel_z_pvt,ddm_utc);

rx_vel_xyz = [rx_vel_x,rx_vel_y,rx_vel_z];

rx_roll_deg = interp1(pvt_utc,rx_roll_deg_pvt,ddm_utc);
rx_pitch_deg = interp1(pvt_utc,rx_pitch_deg_pvt,ddm_utc);
rx_yaw_deg = interp1(pvt_utc,rx_yaw_deg_pvt,ddm_utc);

rx_attitude = [rx_roll_deg,rx_pitch_deg,rx_yaw_deg];

rx_clk_bias_m = interp1(pvt_utc,rx_clk_bias_m_pvt,ddm_utc);
rx_clk_drift_mps = interp1(pvt_utc,rx_clk_drift_mps_pvt,ddm_utc);

rx_clk = [rx_clk_bias_m,rx_clk_drift_mps];
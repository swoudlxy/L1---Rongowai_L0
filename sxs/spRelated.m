% this function computes the sp-related variables, including angles in
% various coordinate frames, ranges, EIRP, nadir antenna gain etc
% Inputs:
% 1) tx, rx: tx and rx structures
% 2) sx_pos_xyz: sx ECEF position vector
% 3) SV_PRN_LUT,SV_eirp_LUT: look-up table between SV number and PRN
% Outputs:
% 1) sp_angle_body: sp angle in body frame, az and theta
% 2) sp_angle_enu: sp angle in ENU frame, az and theta
% 3) theta_gps: GPS off boresight angle
% 4) range: tx to sx range, and rx to sx range
% 5) gps_rad: EIRP, tx power

function [sp_angle_body,sp_angle_enu,sp_angle_ant,theta_gps,range,gps_rad] = spRelated(tx,rx,sx_pos_xyz,SV_eirp_LUT)

% sparse structres
tx_pos_xyz = tx.tx_pos_xyz;
tx_vel_xyz = tx.tx_vel_xyz;
sv_num = tx.sv_num;

rx_pos_xyz = rx.rx_pos_xyz;
rx_vel_xyz = rx.rx_vel_xyz;
rx_att = rx.rx_attitude;

% compute angles
[theta_gps,~] = ecef2orf(tx_pos_xyz,tx_vel_xyz,sx_pos_xyz);

% correct for using radians for ac Euler angles
[sp_theta_body,sp_az_body] = ecef2brf3(rx_pos_xyz,rx_vel_xyz,sx_pos_xyz,rx_att);
[sp_theta_enu,sp_az_enu] = ecef2enuf(rx_pos_xyz,sx_pos_xyz);

sp_angle_body = [sp_theta_body,sp_az_body];
sp_angle_enu = [sp_theta_enu,sp_az_enu];

% compute ranges
R_tsx = norm(sx_pos_xyz-tx_pos_xyz);        % range from tx to sx
R_rsx = norm(sx_pos_xyz-rx_pos_xyz);        % range from rx to sx

range = [R_tsx R_rsx];

% compute gps radiation properties
j = SV_eirp_LUT(:,1) == sv_num;             % index of SV number in eirp LUT

gps_pow_dbw = SV_eirp_LUT(j,3);             % gps power in dBw

% coefficients to compute gps antenna gain
a = SV_eirp_LUT(j,4); b = SV_eirp_LUT(j,5); c = SV_eirp_LUT(j,6); d = SV_eirp_LUT(j,7);
e = SV_eirp_LUT(j,8); f = SV_eirp_LUT(j,9);

% fitting antenna gain
gps_gain_dbi = a*theta_gps^5+b*theta_gps^4+c*theta_gps^3+d*theta_gps^2+e*theta_gps+f;

% compute static gps eirp
stat_eirp_dbw = gps_pow_dbw+gps_gain_dbi;   % static eirp in dbw
stat_eirp_watt = 10^(stat_eirp_dbw/10);     % static eirp in linear watts

gps_rad = [gps_pow_dbw,gps_gain_dbi,stat_eirp_watt];

% compute angles in nadir antenna frame and rx gain
sp_theta_ant = sp_theta_body;
sp_az_ant = sp_az_body+180;

if sp_az_ant > 360
    sp_az_ant = sp_az_ant-360;
end

sp_angle_ant = [sp_theta_ant,sp_az_ant];


















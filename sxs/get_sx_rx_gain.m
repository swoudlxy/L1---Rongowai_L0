function sx_rx_gain = get_sx_rx_gain(sp_angle_ant,nadir_pattern)

% define azimuth and elevation angle in the antenna frame
res = 0.1;                                          % resolution in degrees
az_deg = 0:res:360; 
el_deg = 120:-1*res:0;

lhcp_gain_pattern = nadir_pattern.LHCP;
rhcp_gain_pattern = nadir_pattern.RHCP;

sp_theta_ant = sp_angle_ant(1);
sp_az_ant = sp_angle_ant(2);

[~,az_index] = min(abs(sp_az_ant-az_deg));
[~,el_index] = min(abs(sp_theta_ant-el_deg));

lhcp_gain_dbi = lhcp_gain_pattern(az_index,el_index);
rhcp_gain_dbi = rhcp_gain_pattern(az_index,el_index);

sx_rx_gain = [lhcp_gain_dbi,rhcp_gain_dbi];
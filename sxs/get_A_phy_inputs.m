% this function get three main input variable range for the a A_phy LUT

function [rx_alt_bins,inc_angle_bins,az_angle_bins] = get_A_phy_inputs(input_file)

fid = fopen(input_file);

min_rx_alt = fread(fid,1,'uint16');
res_rx_alt = fread(fid,1,'uint16');
num_rx_alt = fread(fid,1,'uint16');

min_inc_angle = fread(fid,1,'uint16');
res_inc_angle = fread(fid,1,'uint16');
num_inc_angle = fread(fid,1,'uint16');

min_az_angle = fread(fid,1,'uint16');
res_az_angle = fread(fid,1,'uint16');
num_az_angle = fread(fid,1,'uint16');

rx_alt_bins = min_rx_alt+((1:1:num_rx_alt)-1)*res_rx_alt;
inc_angle_bins = min_inc_angle+((1:1:num_inc_angle)-1)*res_inc_angle;
az_angle_bins = min_az_angle+((1:1:num_az_angle)-1)*res_az_angle;
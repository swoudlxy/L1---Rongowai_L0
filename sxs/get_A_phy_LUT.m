% this function retrieves A_phy LUT and the input matrices

function [rx_alt_bins,inc_angle_bins,az_angle_bins,A_phy_LUT_all] = get_A_phy_LUT(LUT_file)

fid = fopen(LUT_file);

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

A_phy_LUT_all = zeros(num_rx_alt,num_inc_angle,num_az_angle,7,41)+nan;

for m = 1:num_rx_alt
    for n = 1:num_inc_angle
        for k = 1:num_az_angle

            data = fread(fid,7*41,'uint32');
            data = reshape(data,[7 41]);

            A_phy_LUT_all(m,n,k,:,:) = data;


        end
    end
end

fclose(fid);
% this function retrieves A_phy LUT

function A_phy_LUT = get_A_phy_LUT(LUT_file,num_rx_alt,num_inc_angle,num_az_angle)

fid = fopen(LUT_file,'r');

A_phy_LUT = zeros(num_rx_alt,num_inc_angle,num_az_angle,7,41)+nan;

for m = 1:num_rx_alt
    for n = 1:num_inc_angle
        for k = 1:num_az_angle

            data = fread(fid,7*41,'uint32');
            data = reshape(data,[7 41]);

            A_phy_LUT(m,n,k,:,:) = data;

        end
    end
end

fclose(fid);
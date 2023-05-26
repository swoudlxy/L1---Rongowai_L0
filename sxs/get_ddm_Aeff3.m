
function [A_eff,A_eff_all,A_phy] = get_ddm_Aeff3(tx,rx,sx, ...
    delay_res,doppler_res, ...
    local_dem,phy_ele_size,chi2)

% sparse structures
tx_pos_xyz = tx.tx_pos_xyz;
tx_vel_xyz = tx.tx_vel_xyz;

rx_pos_xyz = rx.rx_pos_xyz;
rx_vel_xyz = rx.rx_vel_xyz;

sx_pos_xyz = sx.sx_pos_xyz;
sx_pos_lla = ecef2lla(sx_pos_xyz);

sx_delay_bin = sx.sx_delay_bin+1;               % need to fix all 0-indexed bin to 1-indexed
sx_doppler_bin = sx.sx_doppler_bin+1;

% sparse local_dem structure
local_lat = local_dem.lat;
local_lon = local_dem.lon;
local_ele = local_dem.ele;

num_grids = length(local_lat);

% delay & Doppler at SP
[delay_chips_sx,doppler_hz_sx] = deldop(tx_pos_xyz,rx_pos_xyz, ...
    tx_vel_xyz,rx_vel_xyz,sx_pos_xyz);

delay_chips_local = zeros(num_grids)+nan;
doppler_hz_local = zeros(num_grids)+nan;

for m = 1:num_grids
    for n = 1:num_grids

        p_lla1 = [local_lat(m) local_lon(n) local_ele(m,n)];
        p_pos_xyz1 = lla2ecef(p_lla1);

        [delay_p1,doppler_p1,~] = deldop(tx_pos_xyz,rx_pos_xyz, ...
            tx_vel_xyz,rx_vel_xyz,p_pos_xyz1);

        delay_chips_local(m,n) = delay_p1-delay_chips_sx;
        doppler_hz_local(m,n) = doppler_p1-doppler_hz_sx;

    end
end

% interpolate to a fine resolution
sample_rate = 5;                                % dense rate
num_grids_fine = (num_grids-1)*sample_rate+1;   % densed number of grids

local_lat_fine = linspace(local_lat(1),local_lat(end),num_grids_fine);
local_lon_fine = linspace(local_lon(1),local_lon(end),num_grids_fine);

delay_chips = interp2(local_lon,local_lat',delay_chips_local, ...
    local_lon_fine,local_lat_fine','spline');
doppler_hz = interp2(local_lon,local_lat',doppler_hz_local, ...
    local_lon_fine,local_lat_fine','spline');

% get physical size
sx_pos_lat = sx_pos_lla(1);
[~,idx_lat] = min(abs(phy_ele_size(:,1)-sx_pos_lat));
dA = phy_ele_size(idx_lat-floor(num_grids_fine/2):idx_lat+floor(num_grids_fine/2),2);

dA = repmat(dA,1,num_grids_fine);

% construct physical size DDM
A_phy = zeros(5,40);

% bin to physical size DDM
for m = 1:num_grids_fine
    for n = 1:num_grids_fine

        delay_bin_idx1 = floor(-1*delay_chips(m,n)/delay_res+sx_delay_bin);
        doppler_bin_idx1 = floor(doppler_hz(m,n)/doppler_res+sx_doppler_bin);

        if (delay_bin_idx1>=1) && (delay_bin_idx1<=40) && ...
                (doppler_bin_idx1>=1) && (doppler_bin_idx1<=5)

            temp = A_phy(doppler_bin_idx1,delay_bin_idx1);
            temp = temp+dA(m,n);
            A_phy(doppler_bin_idx1,delay_bin_idx1) = temp;

        end

    end
end

% convolution to A_eff
A_eff_all = conv2(A_phy,chi2);
A_eff = A_eff_all(3:7,20:59);          % cut suitable size for A_eff
%A_eff = A_eff_all(3:7,1:40);

% this function computes the effective scattering area at the given surface
% Inputs:
% 1) tx, rx: tx and rx structures
% 2) sx_pos_xyz: ecef position of specular points
% 3) ddm: ddm structure
% 4) local_dem: local region centred at sx
% 5) T_coh: coherent integration duration
% Output:
% 1) A_eff: effective scattering area
% 2) sp_delay_bin,sp_doppler_bin: floating specular bin

function [A_eff,A_eff_all,A_phy] = get_ddm_Aeff(tx,rx,sx,local_dem,phy_ele_size,chi2)

delay_res = 1023000/8192000;            % corrected delay resolution
doppler_res = 250;                      % corrected doppler resolution                        

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

% get coarsen local_dem
sample_rate = 30;

lat_coarse = local_lat(1:sample_rate:end);
lon_coarse = local_lon(1:sample_rate:end);
ele_coarse = local_ele(1:sample_rate:end,1:sample_rate:end);

num_grids_coarse = length(lat_coarse);

% get delay-doppler map over the surface
delay_coarse = zeros(num_grids_coarse);        
doppler_coarse = zeros(num_grids_coarse);

[delay_chips_sx,doppler_Hz_sx] = deldop(tx_pos_xyz,rx_pos_xyz, ...
    tx_vel_xyz,rx_vel_xyz,sx_pos_xyz);

for m = 1:num_grids_coarse
    for n = 1:num_grids_coarse
                    
        p_pos_lla1 = [lat_coarse(m) lon_coarse(n) ele_coarse(m,n)];
        p_pos_xyz1 = lla2ecef(p_pos_lla1);
        
        [delay_p1,doppler_p1,~] = deldop(tx_pos_xyz,rx_pos_xyz, ...
            tx_vel_xyz,rx_vel_xyz,p_pos_xyz1);

        delay_coarse(m,n) = delay_p1-delay_chips_sx;
        doppler_coarse(m,n) = doppler_p1-doppler_Hz_sx;
                                        
    end
end

% interpolate to 30-m resolution
delay_chips = interp2(lon_coarse,lat_coarse',delay_coarse, ...
    local_lon,local_lat','spline');
doppler_Hz = interp2(lon_coarse,lat_coarse',doppler_coarse, ...
    local_lon,local_lat','spline');

% get physical size
sx_pos_lat = sx_pos_lla(1);
[~,idx_lat] = min(abs(phy_ele_size(:,1)-sx_pos_lat));
dA = phy_ele_size(idx_lat-floor(num_grids/2):idx_lat+floor(num_grids/2),2);

dA = repmat(dA,1,num_grids);

% construct physical size DDM
A_phy = zeros(5,40);
     
% bin to physical size DDM
for m = 1:num_grids
    for n = 1:num_grids

        delay_bin_idx1 = floor(-1*delay_chips(m,n)/delay_res+sx_delay_bin);
        doppler_bin_idx1 = floor(doppler_Hz(m,n)/doppler_res+sx_doppler_bin);

        if (delay_bin_idx1>=1) && (delay_bin_idx1<=40) && ...
                (doppler_bin_idx1>=1) && (doppler_bin_idx1<=5)

            temp = A_phy(doppler_bin_idx1,delay_bin_idx1);
            temp = temp+dA(m,n);
            A_phy(doppler_bin_idx1,delay_bin_idx1) = temp;

        end

    end
end

% convolution to A_eff
A_eff1 = conv2(A_phy,chi2);
A_eff = A_eff1(3:7,20:59);          % cut suitable size for A_eff
%A_eff = A_eff1(3:7,1:40);
A_eff_all = A_eff1;
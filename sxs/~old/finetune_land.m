% This function fine-tunes the SP location over land DEMs
% multiple peak processing disabled
% Inputs
% 1)tx - tx structure, including tx position and velocity
% 2)rx - rx structure, including rx position, velocity, and rx clock drifts
% 3)sx_lla_coarse - sx coordinate on WGS84
% 4)ddm - ddm structure
% 5)dem - land dem, including lat, lon and elevation
% 6)L_m, res_m - side length of searching area, resolution of DEM
% Outputs
% 1) sx_xyz_final - fine-tuned sp coordinates
% 2) theta_i - local incidence angle

%function [sx_xyz_final,theta_i] = finetune_land(tx,rx,sx_lla_coarse,L_m,res_m,ddm,dem1,dem2)
clear
clc

% load SRTM_30 DEM
dem_path = '../dat/dem/';
dem_file1 = 'nzsrtm_30_part1_v1.dat';
dem_file2 = 'nzsrtm_30_part2_v1.dat';

dem1 = get_dem([dem_path dem_file1]);
dem2 = get_dem([dem_path dem_file2]);

% load DTU10 model
dtu_path = '../dat/dtu/';
dtu_filename = 'dtu10_v1.dat';

dtu10 = get_dtu10([dtu_path dtu_filename]);

% load distance to coast mask
landmask_path = '../dat/cst/';
landmask_filename = 'dist_to_coast_nz_v1.dat';

dist_to_coast_nz = dist_to_coast_mask([landmask_path landmask_filename]);

%%
clc

load('../exp/tx1.mat');
load('../exp/rx1.mat');
load('../exp/ddm1.mat');

tx = tx1;   rx = rx1;   ddm = ddm1;

clear tx1 rx1 ddm1

sx_lla_coarse = [-43.2688167919431	172.533614183098	0];

L_m = 6030; res_m = 30;

%% code starts from here
clc

% sparse tx structure
tx_pos_xyz = tx.tx_pos_xyz;
tx_vel_xyz = tx.tx_vel_xyz;

% sparse rx structure
rx_pos_xyz = rx.rx_pos_xyz;
rx_vel_xyz = rx.rx_vel_xyz;
rx_clk_drift = rx.rx_clk_drift;

% sparse ddm structure
% rotate raw_count matrix so row is delay and col is doppler
raw_counts1 = ddm.raw_counts;
raw_counts = raw_counts1';                      

delay_dir_chips = ddm.delay_dir_chips;

delay_res = ddm.delay_bin_res;
doppler_res = ddm.doppler_bin_res;

delay_center_bin = ddm.delay_center_bin;
doppler_center_bin = ddm.doppler_center_bin;

delay_center_chips = ddm.delay_center_chips;
doppler_center_Hz = ddm.doppler_center_Hz;

% retrieve a local elevation of the searching region centred at sp_coarse
search_region = get_local_dem(sx_lla_coarse,L_m,res_m,dem1,dem2);
search_lat = search_region.lat;
search_lon = search_region.lon;
search_ele = search_region.ele;

% dimension of the searching area
L = length(search_lat);
M = L-2;    N = L-2;

% contruct <lat/lon> matrix of the search region
[search_lat_matrix,search_lon_matrix] = meshgrid(search_lat(2:end-1),search_lon(2:end-1));

% delay and Doppler bin of the peak power at raw counts
peak_raw_counts = max(max(raw_counts));            % peak power
[delay_peak_bin,doppler_peak_bin] = find(raw_counts == peak_raw_counts,1,'first');

I = length(delay_peak_bin);                         % number of peaks

% compute delay, Doppler values of the searching area
% compute snell angle difference of the searching area

% initialise variables
add_delay_chips = zeros(L-2);                       % additonal range delay in chips
doppler_Hz = zeros(L-2);                            % absolute Doppler frequency in Hz
theta_i_all = zeros(L-2);                           % incidence angle over the searching area

d_delay_chip = zeros(L-2);                          % difference between observed and estimated delay
d_doppler_Hz = zeros(L-2);                          % difference between observed and estimated Doppler
d_snell_deg = zeros(L-2);                           % difference of Snell angle

for m = 2:L-1
    for n = 2:L-1
        
        % coordinate and elevation of the pixel
        % 3*3 matrix is used to compute local incidence angle
        lat1 = search_lat(m-1:m+1);         lat_c = lat1(2);
        lon1 = search_lon(n-1:n+1);         lon_c = lon1(2);
        ele1 = search_ele(m-1:m+1,n-1:n+1); ele_c = ele1(2,2);

        % delay, Doppler and additional phase delay of the pixel
        [~,doppler1,add_delay_chips1] = deldop(tx_pos_xyz,rx_pos_xyz,tx_vel_xyz,rx_vel_xyz, ...
            lat_c,lon_c,ele_c);

        % criteria 3 - Snell reflection law
        % compute theta and phi in degrees over the searching area
        [theta_i1,theta_s1,phi_i1,phi_s1] = angles(lat1,lon1,ele1,tx_pos_xyz,rx_pos_xyz);

        % difference of theta_i and theta_s, and phi_i and phi_s in degrees
        d_theta1 = theta_i1-theta_s1;
        d_phi2 = sind(phi_s1-(phi_i1+180))/cosd(phi_s1-(phi_i1+180));
        d_phi1 = atand(d_phi2);

        % angle difference
        d_snell1 = abs(d_theta1)+abs(d_phi1);

        add_delay_chips(m-1,n-1) = add_delay_chips1;
        doppler_Hz(m-1,n-1) = doppler1;
        d_snell_deg(m-1,n-1) = d_snell1;
        theta_i_all(m-1,n-1) = theta_i1;

    end
end

% loop over the searching area to determine the SP's delay
% and doppler against criteria 1 - minimal delay difference, and 
% 2 - minimal doppler difference

% identify NGRx-observed absolute delay and doppler values for each DDM
% sample (delay_max, and doppler_max)
delay_max_chips = delay_center_chips+delay_res*(delay_peak_bin-delay_center_bin);
doppler_max_Hz = doppler_center_Hz+doppler_res*(doppler_peak_bin-doppler_center_bin);

% loop over the surface
for m = 1:M
    for n = 1:N

        % additional code phase delay and doppler for pixel (m,n)
        add_delay_chips1 = add_delay_chips(m,n);
        doppler1 = doppler_Hz(m,n);

        % criteria 1 - minimal delay/chip difference
        delay_ref2 = delay_dir_chips-add_delay_chips1;  % reference code phase delay for pixel (m,n)
        delay_ref1 = delay_correction(delay_ref2);
            
        % delay difference
        d_delay1 = delay_max_chips-delay_ref1;
        
        % criteria 2 - minimal Doppler difference   
        c = 299792458;
        doppler_clk = rx_clk_drift/c;                   % doppler due to rx drifts
        doppler1_adj = doppler1+doppler_clk;            % adjusted Doppler
            
        % doppler difference
        d_doppler1 = doppler_max_Hz-doppler1_adj;

        % save to matrices
        d_delay_chip(m,n) = d_delay1;
        d_doppler_Hz(m,n) = abs(d_doppler1);
        
    end
end
        
% find pixels meet all three criteria
% these three criteria are defined based on CYGNSS data, need to be
% refined for Rongowai
index = find(abs(d_delay_chip) <= 2.5 & d_doppler_Hz <= 200 & d_snell_deg <= 2);

lat_final = search_lat_matrix(index);
lon_final = search_lon_matrix(index);

%%
[~,lat_final_index] = min(abs(search_lat-lat_final)); 
[~,lon_final_index] = min(abs(search_lon-lon_final));
ele_final = search_ele(lat_final_index,lon_final_index);

    sx_lla_final = [lat_final,lon_final,ele_final];
    sx_xyz_final = lla2ecef([lat_final,lon_final,ele_final]);

    theta_i_final = theta_i_all(lat_final_index,lon_final_index);
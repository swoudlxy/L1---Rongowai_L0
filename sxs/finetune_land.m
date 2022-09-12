% This function fine-tunes the SP location over land DEMs
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

function [sx_xyz_final,theta_i] = finetune_land(tx,rx,sx_lla_coarse,ddm,dem,L_m,res_m)

% sparse tx structure
tx_pos_xyz = tx.tx_pos_xyz;
tx_vel_xyz = tx.tx_vel_xyz;

% sparse rx structure
rx_pos_xyz = rx.rx_pos_xyz;
rx_vel_xyz = rx.rx_vel_xyz;
rx_clock_drift = rx.rx_clock_drift_mps;

% sparse dem structure
lat_dem = dem.lat;
lon_dem = dem.lon;
ele_dem = dem.ele;

% sparse ddm structure
% rotate raw_count matrix so row - delay and col - doppler according to CYGNSS
% need to double check with real Rongowai data
raw_counts1 = ddm.raw_counts;
raw_counts = raw_counts1';                      

delay_dir_chips = ddm.delay_dir_chips;

delay_res = ddm.delay_bin_res;
doppler_res = ddm.doppler_bin_res;

delay_center_bin = ddm.delay_center_bin;
doppler_center_bin = ddm.doppler_center_bin;

delay_center_chips = ddm.delay_center_chips;
doppler_center_Hz = ddm.doppler_center_Hz;

% define L - number of grids (single sidelength) of the search area
L = ceil(L_m/res_m); centre_L = ceil(L/2);
M = L-2; N = L-2;

% retrieve the lat/lon at sp_coarse
lat0 = sx_lla_coarse(1); [~,lat0_row] = min(abs(lat0-lat_dem));
lon0 = sx_lla_coarse(2); [~,lon0_col] = min(abs(lon0-lon_dem));

% retrieve elevation of the searching region centred at sp_coarse
search_lat = lat_dem(lat0_row-centre_L:lat0_row+centre_L);
search_lon = lon_dem(lon0_col-centre_L:lon0_col+centre_L);
search_ele = ele_dem(lat0_row-centre_L:lat0_row+centre_L,lon0_col-centre_L:lon0_col+centre_L);

% contruct <lat/lon> matrix of the search region
[search_lat_matrix,search_lon_matrix] = meshgrid(search_lat(2:end-1),search_lon(2:end-1));

% delay and Doppler bin of the peak power at raw counts
% a ddm may have single or multiple peaks depending on the Earth surface
% the code can process both scenarios, but can start from simpler signle peak processing 

% single peak - treat all raw-count ddms with only single peak bin
% if two bins have the same peak power value, take the first bin
% normally one peak corresponding to one specular reflection
% SP solver first implement the single-peak SP solution and onto multiple
% peak solution after getting more data and understanding the frequency of
% flying over mountainous areas
peak_raw_counts = max(max(raw_counts));            % peak power
[delay_peak_bin,doppler_peak_bin] = find(raw_counts == peak_raw_counts,1,'first');

%{
% multiple peaks - multi-peak detection is performed for multi-sp detection
% identify peak power bins 98% -100%  of the maximal power
% the range needs to be identified after experimenting with more real data
peak_counts_max = max(raw_counts,[],'all');         % peak power - upper boundary
peak_counts_min = 0.98*peak_counts_max;             % peak power - lower boundary

[delay_peak_bin1,doppler_peak_bin1] = find(raw_counts<=peak_counts_max & raw_counts>=peak_counts_min);
I = length(delay_peak_bin1);                        % number of peaks

delay_peak_bin2 = zeros(I,1);
doppler_peak_bin2 = zeros(I,1);

% detect amount of peaks only if two peak bins are not neighbouring bins
for i = 1:I-1

    delay_bin1 = delay_peak_bin1(i);        delay_bin2 = delay_peak_bin1(i+1);
    doppler_bin1 = doppler_peak_bin1(i);    doppler_bin2 = doppler_peak_bin1(i+1);

    if (delay_bin1==delay_bin2) && (doppler_bin1==doppler_bin2-1)

        peak1 = raw_counts(delay_bin1,doppler_bin1);
        peak2 = raw_counts(delay_bin2,doppler_bin2);

        if peak1>=peak2
            delay_peak_bin2(i) = delay_bin1; doppler_peak_bin2(i) = doppler_bin1;

        else
            delay_peak_bin2(i) = delay_bin2; doppler_peak_bin2(i) = doppler_bin2;

        end
        
    elseif (delay_bin1==delay_bin2-1) && (doppler_bin1==doppler_bin2)

        peak1 = raw_counts(delay_bin1,doppler_bin1);
        peak2 = raw_counts(delay_bin2,doppler_bin2);

        if peak1>=peak2
            delay_peak_bin2(i) = delay_bin1; doppler_peak_bin2(i) = doppler_bin1;

        else
            delay_peak_bin2(i) = delay_bin2; doppler_peak_bin2(i) = doppler_bin2;

        end

    else
        delay_peak_bin2(i) = delay_bin1;    delay_peak_bin2(i+1) = delay_bin2;
        doppler_peak_bin2(i) = doppler_bin1;  doppler_peak_bin2(i+1) = doppler_bin2;

    end

end

is_zero = ~delay_peak_bin2;
delay_peak_bin = delay_peak_bin2(~is_zero); 
doppler_peak_bin = doppler_peak_bin2(~is_zero);
%}

I = length(delay_peak_bin);                         % number of peaks

% compute delay, Doppler values of the searching area
% compute snell angle difference of the searching area

% initialise variables
add_delay_chip = zeros(L-2);                        % additonal range delay in chips
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

        add_delay_chip(m-1,n-1) = add_delay_chips1;
        doppler_Hz(m-1,n-1) = doppler1;
        d_snell_deg(m-1,n-1) = d_snell1;
        theta_i_all(m-1,n-1) = theta_i1;

    end
end

sx_lla_final = zeros(I,3);
sx_xyz_final = zeros(I,3);
theta_i = zeros(I,1);

% for each peak, loop over the searching area to determine the SP's delay
% and doppler against criteria 1 - minimal delay difference, and 
% 2 - minimal doppler difference
for i = 1:I

    delay_peak_bin = delay_peak_bin(i);
    doppler_peak_bin = doppler_peak_bin(i);

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
            delay_ref1 = delay_dir_chips-add_delay_chips1;  % reference code phase delay for pixel (m,n)

            % adjust estimated phase delay within the range: [0 1023]
            while delay_ref1 < 0
                delay_ref1 = delay_ref1+1023;
            end
            
            % delay difference
            d_delay1 = delay_max_chips-delay_ref1;
        
            % criteria 2 - minimal Doppler difference   
            doppler_clk = rx_clock_drift/c;                 % doppler due to rx drifts
            doppler1_adj = doppler1+doppler_clk;            % adjusted Doppler
            
            % doppler difference
            d_doppler1 = doppler_max_Hz-doppler1_adj;
        
        end
    end
        
    % save to matrices
    d_delay_chip(m,n) = d_delay1;
    d_doppler_Hz(m,n) = abs(d_doppler1);

    % find pixels meet all three criteria
    % these three criteria are defined based on CYGNSS data, need to be
    % refined for Rongowai
    index = find(abs(d_add_tau) <= 2.5 & d_doppler_Hz <= 200 & d_snell_deg <= 2);

    lat_final1 = search_lat_matrix(index);
    lon_final1 = search_lon_matrix(index);

    % define multiple SP locations still need to be refined
    % here treat only one peak DDM - i.e, one SP
    if I == 1
        lat_final = mean(lat_final1);
        lon_final = mean(lon_final1);
    end

    [~,lat_final_index] = min(abs(search_lat-lat_final)); 
    [~,lon_final_index] = min(abs(search_lon-lon_final));
    ele_final = search_ele(lat_final_index,lon_final_index);

    sx_lla_final(i,:) = [lat_final,lon_final,ele_final];
    sx_xyz_final(i,:) = lla2ecef([lat_final,lon_final,ele_final]);

    theta_i_final = theta_i_all(lat_final_index,lon_final_index);
    theta_i(i) = theta_i_final;

end
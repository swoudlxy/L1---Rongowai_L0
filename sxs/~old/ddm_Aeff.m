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

function [A_eff,A_eff_all] = ddm_Aeff(tx,rx,sx,ddm,local_dem,T_coh)

chip_rate = 1.023e6;
tau_c = 1/chip_rate;

% sparse structures
tx_pos_xyz = tx.tx_pos_xyz;
tx_vel_xyz = tx.tx_vel_xyz;

rx_pos_xyz = rx.rx_pos_xyz;
rx_vel_xyz = rx.rx_vel_xyz;

sx_pos_xyz = sx.sx_pos_xyz;
sx_delay_bin_float = sx.sx_delay_bin+1;
sx_doppler_bin_float = sx.sx_doppler_bin+1;

num_delay_bins = ddm.num_delay_bins;
num_doppler_bins = ddm.num_doppler_bins;

delay_center_bin = ddm.delay_center_bin+1;
doppler_center_bin = ddm.doppler_center_bin+1;

d_delay = ddm.delay_resolution;
d_doppler = ddm.doppler_resolution;

sx_delay_bin = floor(sx_delay_bin_float);
delayOffset_frac = sx_delay_bin_float-sx_delay_bin;

sx_doppler_bin = floor(sx_doppler_bin_float);
dopplerOffset_frac = sx_doppler_bin_float-sx_doppler_bin;

% sparse local_dem structure
lat_local = local_dem.lat;
lon_local = local_dem.lon;
ele_local = local_dem.ele;

num_grid = length(lat_local);

% compute physical size of each area element
dA = zeros(num_grid-2);

for m = 2:num_grid-1
    for n = 2:num_grid-1
        
        point1 = [lat_local(m),lon_local(n)];
        
        % coordinates of four points around point1
        point2 = [lat_local(m-1),lon_local(n)];
        point3 = [lat_local(m),lon_local(n-1)];
        point4 = [lat_local(m+1),lon_local(n)];
        point5 = [lat_local(m),lon_local(n+1)];

        % side length in km
        sidelength1 = m_lldist([point1(2) point2(2)],[point1(1) point2(1)]) ... 
            + m_lldist([point1(2) point4(2)],[point1(1) point4(1)]);
        sidelength2 = m_lldist([point1(2) point3(2)],[point1(1) point3(1)]) ... 
            + m_lldist([point1(2) point5(2)],[point1(1) point5(1)]);

        dA(m-1,n-1) = sidelength1*sidelength2*1e6/4;    % size in meter^2

    end
end

% delay and doppler at sx
[delay_chips_sx,doppler_Hz_sx,~] = deldop(tx_pos_xyz,rx_pos_xyz, ...
        tx_vel_xyz,rx_vel_xyz,sx_pos_xyz);

% initalise physical area DDM
DDM_A = zeros(num_doppler_bins,num_delay_bins);

% compute delay and Doppler values over the discretised surface and map
% to physical scattering area
num_grid = num_grid-2;

delay_map = zeros(num_grid);
doppler_map = zeros(num_grid);

for m = 1:num_grid
    for n = 1:num_grid

        % coordinates of the grid
        grid_lat = lat_local(m);
        grid_lon = lon_local(n);
        grid_ele = ele_local(m,n);        

        grid_pos_xyz = lla2ecef([grid_lat grid_lon grid_ele]);

        % absolute delay and doppler values
        [abs_delay_chips,abs_doppler_Hz,~] = deldop(tx_pos_xyz,rx_pos_xyz, ...
              tx_vel_xyz,rx_vel_xyz,grid_pos_xyz);

        % delay and doppler values relative to SX
        grid_delay_chips = abs_delay_chips-delay_chips_sx;
        grid_doppler_Hz = abs_doppler_Hz-doppler_Hz_sx;

        delay_map(m,n) = grid_delay_chips;
        doppler_map(m,n) = grid_doppler_Hz;

        % delay and doppler bin
        delay_bin = round(grid_delay_chips/d_delay+sx_delay_bin+delayOffset_frac);
        doppler_bin = round(grid_doppler_Hz/d_doppler+sx_doppler_bin+dopplerOffset_frac);

        % physical scattering area mapping to DDM
        if (delay_bin<=num_delay_bins) && (delay_bin>0) && ...
                (doppler_bin<=num_doppler_bins) && (doppler_bin>0)

            temp = DDM_A(doppler_bin,delay_bin);
            dA_temp = dA(m,n);
            temp = temp+dA_temp;
            
            DDM_A(doppler_bin,delay_bin) = temp;
            
        end

    end
end

% initialise AF chi
chi = zeros(num_doppler_bins,num_delay_bins);

% compute ambiguity Function (AF, i.e., chi)
for i = 1:num_delay_bins
    for j = 1:num_doppler_bins

        dtau = (i-delay_center_bin)*d_delay*tau_c;      % dtau in second
        dfreq = (j-doppler_center_bin)*d_doppler;       % dfreq in Hz

        % compute complex AF value at each delay-doppler bin
        chi(j,i) = amb_fun(dtau,dfreq,tau_c,T_coh);
    
    end
end

chi_mag = abs(chi);             % magnitude
chi2 = chi_mag.*chi_mag;        % chi_square

temp1 = conv2(chi2,DDM_A);      % 2D convlution
A_eff = temp1(3:7,21:60);       % crop to proper size
A_eff_all = temp1;
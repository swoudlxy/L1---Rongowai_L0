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

function [A_eff,sp_delay_bin_float,sp_doppler_bin_float] = ddm_Aeff(tx,rx,sx_pos_xyz,ddm,local_dem,T_coh)

c = 299792458;                  % light speed metre per second
fc = 1575.42e6;                 % L1 carrier frequency in Hz
lambda = c/fc;                  % wavelength

chip_rate = 1.023e6;            % L1 GPS chip-per-second, code modulation frequency
tau_c = 1/chip_rate;            % C/A code chiping period

% sparse structures
tx_pos_xyz = tx.tx_pos_xyz;
tx_vel_xyz = tx.tx_vel_xyz;

rx_pos_xyz = rx.rx_pos_xyz;
rx_vel_xyz = rx.rx_vel_xyz;
rx_clock_drift_mps = rx.rx_clock_drift_mps;

delay_dir_chips = ddm.delay_dir_chips;

d_delay = ddm.delay_bin_res; 
d_doppler = ddm.doppler_bin_res;

delay_center_bin = ddm.delay_center_bin;
doppler_center_bin = ddm.doppler_center_bin;

delay_center_chips = ddm.delay_center_chips;
doppler_center_Hz = ddm.doppler_center_Hz;

num_delay_bin = ddm.num_delay_bin;
num_doppler_bin = ddm.num_doppler_bin;

lat_local = local_dem.lat_local;
lon_local = local_dem.lon_local;
ele_local = local_dem.ele_local;
dA = local_dem.dA;

% floating specular bin
% additional path phase delay according to the computed SP for 
% multiple SPs, using the 1st SP as the primary SP
R_tsx = norm(sx_pos_xyz(1,:)-tx_pos_xyz);
R_rsx = norm(sx_pos_xyz(1,:)-rx_pos_xyz);
R_trx = norm(rx_pos_xyz-tx_pos_xyz);

% ground-estimated additional range to SP
add_path_delay_chips = meter2chips((R_tsx+R_rsx)-R_trx);
sp_bin_delay_chips = delay_dir_chips+add_path_delay_chips;

while sp_bin_delay_chips > 1023
    sp_bin_delay_chips = mod(sp_bin_delay_chips,1023);
end

% sp delay relative to ddm reference and floating delay bin
d_sp_delay_chips = delay_center_chips-sp_bin_delay_chips;
sp_delay_bin_float = delay_center_bin-d_sp_delay_chips/d_delay;

sp_delay_bin = floor(sp_delay_bin_float);
delayOffset_frac = mod(sp_delay_bin_float,floor(sp_delay_bin_float));

% absolute doppler at SP
tsx_vector = sx_pos_xyz(1,:)-tx_pos_xyz; 
tsx_unit = tsx_vector/norm(tsx_vector);

srx_vector = rx_pos_xyz-sx_pos_xyz(1,:); 
srx_unit = srx_vector/norm(srx_vector);

term1 = dot(tx_vel_xyz,tsx_unit);
term2 = dot(rx_vel_xyz,srx_unit);
term3 = rx_clock_drift_mps/c;

doppler_temp = (term1-term2)/lambda;
sp_bin_doppler_Hz = doppler_temp+term3;

% sp doppler relative to ddm reference and floating doppler bin
d_sp_doppler_Hz = doppler_center_Hz-sp_bin_doppler_Hz;
sp_doppler_bin_float = doppler_center_bin-d_sp_doppler_Hz/d_doppler;

sp_doppler_bin = floor(sp_doppler_bin_float);
dopplerOffset_frac = mod(sp_doppler_bin_float,floor(sp_doppler_bin_float));

% delay and doppler at sx
sx_lla = ecef2lla(sx_pos_xyz);
[delay_chips_sx,doppler_Hz_sx,~] = deldop(tx_pos_xyz,rx_pos_xyz, ...
        tx_vel_xyz,rx_vel_xyz,sx_lla(1),sx_lla(2),sx_lla(3));

% initalise physical area
DDM_A = zeros(num_doppler_bin,num_delay_bin);

% compute delay and Doppler values over the discretised surface and map
% to physical scattering area
for m = 1:num_grid
    for n = 1:num_grid

        % coordinates of the grid
        grid_lat = lat_local(m);
        grid_lon = lon_local(n);
        grid_ele = ele_local(m,n);

        % absolute delay and doppler values
        [abs_delay_chips,abs_doppler_Hz,~] = deldop(tx_pos_xyz,rx_pos_xyz, ...
              tx_vel_xyz,rx_vel_xyz,grid_lat,grid_lon,grid_ele);

        % delay and doppler values relative to SX
        grid_delay_chips = abs_delay_chips-delay_chips_sx;
        grid_doppler_Hz = abs_doppler_Hz-doppler_Hz_sx;

        % delay and doppler bin
        delay_bin = round(grid_delay_chips/d_delay+sp_delay_bin+delayOffset_frac);
        doppler_bin = round(grid_doppler_Hz/d_doppler+sp_doppler_bin+dopplerOffset_frac);

        % physical scattering area mapping to DDM
        if (delay_bin<=num_delay_bin) && (delay_bin>0) && ...
                (doppler_bin<=num_doppler_bin) && (doppler_bin>0)

            temp = DDM_A(doppler_bin,delay_bin);
            dA_temp = dA(m,n);
            temp = temp+dA_temp;
            
            DDM_A(doppler_bin,delay_bin) = temp;
            
        end

    end
end

% initialise chi
chi = zeros(num_delay_bin,num_doppler_bin);

% compute ambiguity Function (AF, i.e., chi)
for i = 1:num_delay_bin
    for j = 1:num_doppler_bin

        dtau = (i-delay_center_bin)*d_delay*tau_c;      % dtau in second
        dfreq = (j-doppler_center_bin)*d_doppler;       % dfreq in Hz

        % compute complex AF value at each delay-doppler bin
        chi(j,i) = amb_fun(dtau,dfreq,tau_c,T_coh);
    
    end
end

chi_mag = abs(chi);             % magnitude
chi2 = chi_mag.*chi_mag;        % chi_square

temp1 = conv2(chi2,DDM_A);      % 2D convlution
A_eff = temp1(4:8,21:60);       % crop to proper size
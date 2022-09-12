% This function computes absolute delay and doppler values for a given
% pixel whose coordinate is <lat,lon,ele>
% The ECEF position and velocity vectors of tx and rx are also required
% Inputs:
% 1) tx_xyz, rx_xyz: ecef position of tx, rx
% 2) tx_vel, rx_vel: ecef velocity of tx, rx
% 3) lat/lon/elevation of the pixel under computation
% Outputs:
% 1) delay_chips: delay measured in chips
% 2) doppler_Hz: doppler measured in Hz
% 3) add_delay_chips: additional delay measured in chips

function [delay_chips,doppler_Hz,add_delay_chips] = deldop(tx_pos_xyz,rx_pos_xyz,tx_vel_xyz,rx_vel_xyz,lat,lon,ele)

%  common parameters
c = 299792458;                      %light speed metre per second
fc = 1575.42e6;                     %L1 carrier frequency in Hz
lambda = c/fc;                      %wavelength

% direct delay chips from tx to rx
dir_delay_chips = meter2chips(norm(rx_pos_xyz-tx_pos_xyz));

p_xyz = lla2ecef([lat lon ele]);    %convert to ECEF coordinate

V_tp = p_xyz-tx_pos_xyz; R_tp = norm(V_tp); V_tp_unit = V_tp/R_tp;
V_pr = rx_pos_xyz-p_xyz; R_pr = norm(V_pr); V_pr_unit = V_pr/R_pr;

% absolute delay in chips
delay_chips = meter2chips((R_tp+R_pr));     

% additional delay (difference between reflected and direct path) relative to direct path in chips
add_delay_chips = delay_chips-dir_delay_chips;

% absolute Doppler frequency in Hz
term1 = dot(tx_vel_xyz,V_tp_unit);
term2 = dot(rx_vel_xyz,V_pr_unit);

doppler_Hz = (term1-term2)/lambda;  %Doppler in Hz
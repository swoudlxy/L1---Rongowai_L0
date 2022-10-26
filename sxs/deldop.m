% This function computes absolute delay and doppler values for a given
% pixel whose coordinate is <lat,lon,ele>
% The ECEF position and velocity vectors of tx and rx are also required
% Inputs:
% 1) tx_xyz, rx_xyz: ecef position of tx, rx
% 2) tx_vel, rx_vel: ecef velocity of tx, rx
% 3) p_xyz of the pixel under computation
% Outputs:
% 1) delay_chips: delay measured in chips
% 2) doppler_Hz: doppler measured in Hz
% 3) add_delay_chips: additional delay measured in chips

function [delay_chips,doppler_hz,add_delay_chips] = deldop(tx_pos_xyz,rx_pos_xyz,tx_vel_xyz,rx_vel_xyz,p_xyz)

%  common parameters
c = 299792458;                          % light speed metre per second
fc = 1575.42e6;                         % L1 carrier frequency in Hz
lambda = c/fc;                          % wavelength

V_tp = tx_pos_xyz-p_xyz; R_tp = norm(V_tp); V_tp_unit = V_tp/R_tp;
V_rp = rx_pos_xyz-p_xyz; R_rp = norm(V_rp); V_rp_unit = V_rp/R_rp;
V_tr = tx_pos_xyz-rx_pos_xyz; R_tr = norm(V_tr);

delay = R_tp+R_rp;  delay_chips = meter2chips(delay);
add_delay_chips = meter2chips(R_tp+R_rp-R_tr);

% absolute Doppler frequency in Hz
term1 = dot(tx_vel_xyz,V_tp_unit);
term2 = dot(rx_vel_xyz,V_rp_unit);

doppler_hz = -1*(term1+term2)/lambda;   % Doppler in Hz
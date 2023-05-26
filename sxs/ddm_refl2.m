% this function computes the land reflectivity by implementing the xpol
% antenna gain
% inputs
% 1) power_analog: L1a product, DDM power in watt
% 2) eirp_watt: transmitter eirp in watt
% 3) rx_gain_db_i: receiver antenna gain in the direction of SP, in dBi
% 4) R_tsx, R_rsx: tx to sp range and rx to sp range, in meters
% outputs
% 1) copol and xpol reflectivity

function [refl_copol,refl_xpol] = ddm_refl2(power_analog_LHCP,power_analog_RHCP,eirp_watt,rx_gain_db_i,R_tsx,R_rsx)

% define constants
c = 299792458;                              % speed of light, meter per second
freq = 1575.42e6;                           % GPS L1 operating frequency, Hz
lambda = c/freq;                            % wavelength, meter
lambda2 = lambda*lambda;

rx_gain = db2pow(rx_gain_db_i);             % convert antenna gain to linear form

range = R_tsx+R_rsx;

term1 = (4*pi*range)^2;
term2 = eirp_watt*lambda2;
term3 = term1/term2;

term4 = term3*(rx_gain)^(-1);

refl_copol = term4(1,1)*power_analog_LHCP+term4(1,2)*power_analog_RHCP;
refl_xpol  = term4(2,1)*power_analog_LHCP+term4(2,2)*power_analog_RHCP;
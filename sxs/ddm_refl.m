%this function computes the land reflectivity
%inputs
%1)power_analog: L1a product, DDM power in watt
%2)eirp_watt: transmitter eirp in watt
%3)rx_gain_db_i: receiver antenna gain in the direction of SP, in dBi
%4)TSp, RSp: tx to sp range and rx to sp range, in meters
%outputs
%1)reflectivity
%2)reflectivity peak

function [reflectivity,reflectivity_peak] = ddm_refl(power_analog,eirp_watt,rx_gain_db_i,TSp,RSp)

%define constants
c = 299792458;                              %speed of light, meter per second
freq = 1575.42e6;                           %GPS L1 operating frequency, Hz
lambda = c/freq;                            %wavelength, meter
lambda2 = lambda*lambda;

sp_rx_gain_pow = db2pow(rx_gain_db_i);        %convert antenna gain to linear form

range = TSp+RSp;

term1 = (4*pi*range)^2;
term2 = eirp_watt*sp_rx_gain_pow*lambda2;
term3 = term1/term2;

reflectivity = power_analog*term3;
reflectivity_peak = max(reflectivity,[],'all');
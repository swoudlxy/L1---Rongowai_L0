% this version has copol xpol antenna gain implemented
% this function computes bistatic radar cross section (BRCS) according to
% bistatic radar equation based on the inputs as below
% inputs:
% 1) power_analog: L1a product in watts
% 2) eirp_watt, rx_gain_db_i: gps eirp in watts and rx antenna gain in dBi
% 3) TSx, RSx: Tx to Sx and Rx to Sx ranges
% outputs: 
% 1) brcs: bistatic RCS

function [brcs_copol,brcs_xpol] = ddm_brcs2(power_analog_LHCP,power_analog_RHCP,eirp_watt,rx_gain_db_i,TSx,RSx)

% define constants
c = 299792458;                      % light speed, m/s
f = 1575.42e6;                      % GPS L1 band, Hz
lambda = c/f;                       % wavelength, m
lambda2 = lambda*lambda;

% derive BRCS
rx_gain = db2pow(rx_gain_db_i);     % linear rx gain
range = TSx*RSx;

term1 = 4*pi*(4*pi*range)^2;
term2 = eirp_watt*lambda2;
term3 = term1/term2;

term4 = term3*(rx_gain)^(-1);
        
brcs_copol = term4(1,1)*power_analog_LHCP+term4(1,2)*power_analog_RHCP;
brcs_xpol  = term4(2,1)*power_analog_LHCP+term4(2,2)*power_analog_RHCP;
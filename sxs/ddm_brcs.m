% this function computes bistatic radar cross section (BRCS) according to
% bistatic radar equation based on the inputs as below
% inputs:
% 1) power_analog: L1a product in watts
% 2) eirp_watt, rx_gain_db_i: gps eirp in watts and rx antenna gain in dBi
% 3) TSx, RSx: Tx to Sx and Rx to Sx ranges
% outputs: 
% 1) brcs: bistatic RCS

function brcs = ddm_brcs(power_analog,eirp_watt,rx_gain_db_i,TSx,RSx)

% define constants
c = 299792458;                      % light speed, m/s
f = 1575.42e6;                      % GPS L1 band, Hz
lambda = c/f;                       % wavelength, m
lambda2 = lambda*lambda;

rx_gain = db2pow(rx_gain_db_i);     % linear rx gain

term1 = eirp_watt*rx_gain;

term2_1 = TSx*RSx;
term2 = 1/(term2_1*term2_1);

power_factor = lambda2*term1*term2/((4*pi)^3);

brcs = power_analog/power_factor;
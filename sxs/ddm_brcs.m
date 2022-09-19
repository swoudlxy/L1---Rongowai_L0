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

term1 = (4*pi)^3/lambda2;
term2 = (TSx*RSx)^2;
term3 = eirp_watt*rx_gain;

brcs = power_analog*term1*term2/term3;
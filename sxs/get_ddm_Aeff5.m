% this is to construct effective scattering area from LUTs
% this version considers the floating SP Doppler bin

function A_eff = get_ddm_Aeff5(rx_alt,inc_angle,az_angle, ...
    rx_alt_bins,inc_angle_bins,az_angle_bins, ...
    sp_delay_bin,sp_doppler_bin,chi2,A_phy_LUT_all)

% center delay and doppler bin for full DDM size
% note: not to be used elsewhere
center_delay_bin = 40;
center_doppler_bin = 5;

% integer and fractional part of the SP bin
sp_delay_intg = floor(sp_delay_bin);
sp_delay_frac = sp_delay_bin-sp_delay_intg;

sp_doppler_intg = round(sp_doppler_bin);
sp_doppler_frac = sp_doppler_bin-sp_doppler_intg;

% index for the floating Doppler bin
% which LUT should be used for interploation
k = floor((sp_doppler_frac+0.5)/0.1)+1;

if sp_doppler_frac>=0
    A_phy_LUT1 = A_phy_LUT_all(k).A_phy_LUT;
    A_phy_LUT2 = A_phy_LUT_all(k+1).A_phy_LUT;

elseif sp_doppler_frac<0
    A_phy_LUT2 = A_phy_LUT_all(k).A_phy_LUT;
    A_phy_LUT1 = A_phy_LUT_all(k+1).A_phy_LUT;

end

% derive full scattering area from LUT
A_phy_full_1_1 = interp3D(inc_angle_bins,rx_alt_bins,az_angle_bins,A_phy_LUT1, ...
    inc_angle,rx_alt,az_angle);
A_phy_full_2_1 = [zeros(1,41);A_phy_full_1_1;zeros(1,41)];
A_phy_full_1 = [A_phy_full_2_1 zeros(9,39)];

A_phy_full_1_2 = interp3D(inc_angle_bins,rx_alt_bins,az_angle_bins,A_phy_LUT2, ...
    inc_angle,rx_alt,az_angle);
A_phy_full_2_2 = [zeros(1,41);A_phy_full_1_2;zeros(1,41)];
A_phy_full_2 = [A_phy_full_2_2 zeros(9,39)];

A_phy_full = (1-sp_doppler_frac)*A_phy_full_1+sp_doppler_frac*A_phy_full_2;

% shift A_phy_full according to the floating SP delay bin
A_phy_shift = [A_phy_full(:,end) A_phy_full(:,1:end-1)];
A_phy_shift = (1-sp_delay_frac)*A_phy_full+sp_delay_frac*A_phy_shift;

% crop the A_phy_full to Rongowai 5*40 DDM size
delay_shift_bin = center_delay_bin-sp_delay_intg;
doppler_shift_bin = center_doppler_bin-sp_doppler_intg;

A_phy = A_phy_shift(doppler_shift_bin+1:doppler_shift_bin+5, ...
    delay_shift_bin+1:delay_shift_bin+40);

% convolution with WAF to get A_eff
A_eff_all = conv2(A_phy,chi2);
A_eff = A_eff_all(3:7,21:60);


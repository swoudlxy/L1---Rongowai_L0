% this function gets 2D AF

function chi2 = get_chi2(num_delay_bins,num_doppler_bins, ...
    delay_center_bin,doppler_center_bin, ...
    delay_res,doppler_res)

chip_rate = 1.023e6;
tau_c = 1/chip_rate;
T_coh = 1/1000;

%delay_res = 1023000/8192000;    % corrected delay resolution
%doppler_res = 250;              % corrected doppler resolution

chi = zeros(num_doppler_bins,num_delay_bins);

for i = 1:num_delay_bins
    for j = 1:num_doppler_bins

        dtau = (i-delay_center_bin)*delay_res*tau_c;      % dtau in second
        dfreq = (j-doppler_center_bin)*doppler_res;       % dfreq in Hz

        % compute complex AF value at each delay-doppler bin
        chi(j,i) = get_amb_fun(dtau,dfreq,tau_c,T_coh);
    
    end
end

chi_mag = abs(chi);             % magnitude
chi2 = chi_mag.*chi_mag;        % chi_square
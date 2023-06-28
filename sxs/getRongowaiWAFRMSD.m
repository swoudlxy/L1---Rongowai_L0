function outputs = getRongowaiWAFRMSD(ddm,delay_res_chips,rmsd_delay_span_chips,fft_interpolation_factor,power_vs_delay_noise_floor_rows)

s = size(ddm); % [doppler by delay; 5 by 40]
power_vs_delay          = sum(ddm,1);
power_vs_delay          = interpft(power_vs_delay,fft_interpolation_factor*s(2));
[~,delay_peak_index]    = max(power_vs_delay);
downsample_axis         = delay_peak_index-floor(fft_interpolation_factor*s(2)/2):fft_interpolation_factor:delay_peak_index+floor(fft_interpolation_factor*s(2)/2);
downsample_axis(downsample_axis<1)                      = [];
downsample_axis(downsample_axis>length(power_vs_delay)) = [];
power_vs_delay          = power_vs_delay(downsample_axis);
power_vs_delay          = fliplr(power_vs_delay);
noise_floor_counts      = mean(power_vs_delay(1:power_vs_delay_noise_floor_rows));
power_vs_delay          = power_vs_delay-noise_floor_counts;
snr_db                  = 10*log10(max(power_vs_delay)/noise_floor_counts);
power_vs_delay          = power_vs_delay./max(power_vs_delay);

[~,delay_peak_index]                      = find(power_vs_delay==max(power_vs_delay));
delay_axis_chips                          = 1:1:length(power_vs_delay);
delay_axis_chips(1:delay_peak_index-1)    = delay_axis_chips(1:delay_peak_index-1)-delay_axis_chips(delay_peak_index);
delay_axis_chips(delay_peak_index+1:end)  = delay_axis_chips(delay_peak_index+1:end)-delay_axis_chips(delay_peak_index);
delay_axis_chips(delay_peak_index)        = 0;
delay_axis_chips                          = delay_axis_chips.*delay_res_chips;

rwaf                                      = rongowaiWAF(3,delay_peak_index,[s(1) length(delay_axis_chips)],delay_axis_chips);
rwaf                                      = nansum(rwaf,1);
rwaf                                      = rwaf./nanmax(rwaf);

outputs.delay_axis_chips                  = delay_axis_chips;
outputs.waf                               = rwaf;
outputs.power_vs_delay                    = power_vs_delay;
outputs.power_vs_delay_snr_db             = snr_db;

span_over_which_to_compute_rmsd = delay_peak_index-round(rmsd_delay_span_chips/delay_res_chips):1:delay_peak_index+round(rmsd_delay_span_chips/delay_res_chips);
span_over_which_to_compute_rmsd(span_over_which_to_compute_rmsd<1) = [];
span_over_which_to_compute_rmsd(span_over_which_to_compute_rmsd>length(power_vs_delay)) = [];

rwaf                    = rwaf(span_over_which_to_compute_rmsd);
power_vs_delay          = power_vs_delay(span_over_which_to_compute_rmsd);
rmsd                    = sqrt(nansum(power(rwaf-power_vs_delay,2)));
outputs.rmsd            = rmsd;

end


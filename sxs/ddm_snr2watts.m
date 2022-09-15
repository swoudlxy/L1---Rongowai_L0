% this function converts from ddm_snr_db to power in watts
% this is part of L1a calibration
% Inputs:
% 1) ddm_snr_db: DDM SNR measured in dB
% 2) noise_std_ch: noise standard deviation of that channel
% 3) cal_factors_ch: signal and noise calibration factors in dB of that
% channel
% Outputs:
% 1) signal_power_watts: signal power per bin
% 2) noise_power_watts: noise power per bin
% 3) snr_db: ratio between maximal signal power and average noise power

function [signal_power_watts,noise_power_watts,snr_db] = ddm_snr2watts(ddm_snr_db,noise_std_ch,cal_factors_ch)
%cal_factors_ch = cal_factors_ch2;
%noise_std_ch = noise_std_ch2;

ddm_snr = db2pow(ddm_snr_db);

noise_power_14bit_db = mag2db(noise_std_ch);

noise_power_14bit = db2pow(noise_power_14bit_db);
signal_power_14bit = ddm_snr*noise_power_14bit-noise_power_14bit;
%signal_power_14bit_db = pow2db(ddm_snr*noise_power_14bit-noise_power_14bit);

signal_power_watts = signal_power_14bit*db2pow(cal_factors_ch(1));
noise_power_watts = noise_power_14bit*db2pow(cal_factors_ch(2));

%signal_power_dbm = signal_power_14bit_db+cal_factors_ch(1);
%noise_power_dbm = noise_power_14bit_db+cal_factors_ch(2);

%signal_power_watts = db2pow(signal_power_dbm);
%noise_power_watts = db2pow(noise_power_dbm);

signal_power_peak = max(max(signal_power_watts));
noise_avg = mean(noise_power_watts);
snr_db = (signal_power_peak/noise_avg);
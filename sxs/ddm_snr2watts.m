% this function converts from ddm_snr_db to power in watts
% this is part of L1a calibration
% Inputs:
% 1) ddm_snr_db: DDM SNR measured in dB
% 2) noise_std_ch: noise standard deviation of that channel
% 3) cal_factors_ch: signal and noise calibration factors in dB
% Outputs:
% 1) signal_power_watts: signal power per bin
% 2) noise_power_watts: noise power per bin

function [signal_power_watts,noise_power_watts] = ddm_snr2watts(ddm_snr_db,noise_std_ch,cal_factors_ch)

% all computations are implemented in linear power
ddm_snr = db2pow(ddm_snr_db);

noise_power_14bit = noise_std_ch^2;
signal_power_14bit = ddm_snr*noise_power_14bit-noise_power_14bit;

signal_power_watts = signal_power_14bit*db2pow(cal_factors_ch(1));
noise_power_watts = noise_power_14bit*db2pow(cal_factors_ch(2));
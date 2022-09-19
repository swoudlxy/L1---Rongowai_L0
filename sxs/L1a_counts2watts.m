% this function performs Rongowai L1a calibration by converting from raw
% measured counts to power in watts
% Inputs:
% 1) raw_counts: raw measurments in counts
% 2) ANZ_port: signal source
% 3) noise_std: noise standard deviation
% Outputs:
% 1) signal_power_watts: signal power per bin
% 2) noise_power_watts: noise power per bin
% 3) snr: ratio between maximal signal power and average noise power

function [signal_power_watts,noise_power_watts,snr_db] = L1a_counts2watts(raw_counts,ANZ_port,noise_std)

% convert to db
raw_counts_db = mag2db(raw_counts);
%ddm_snr_db = raw_counts_db/2-78.3;          % 78.3 dB is the noise-only DDM counts
ddm_snr_db = raw_counts_db-78.3;

% define calibration factors
signal_cal_factors_db = [-162.93,-163.98,-161.24];
noise_cal_factors_db = [-133.51,-134.45,-132.13];

% perform calibration based on the input ANZ port channel
switch ANZ_port

    case 1
        noise_std_ch1 = noise_std(1);
        cal_factors_ch1 = [signal_cal_factors_db(1),noise_cal_factors_db(1)];

        [signal_power_watts,noise_power_watts,snr_db] = ddm_snr2watts(ddm_snr_db,noise_std_ch1,cal_factors_ch1);
        
    case 2
        noise_std_ch2 = noise_std(2);
        cal_factors_ch2 = [signal_cal_factors_db(2),noise_cal_factors_db(2)];

        [signal_power_watts,noise_power_watts,snr_db] = ddm_snr2watts(ddm_snr_db,noise_std_ch2,cal_factors_ch2);
            
    case 3
        noise_std_ch3 = noise_std(3);
        cal_factors_ch3 = [signal_cal_factors_db(3),noise_cal_factors_db(3)];

        [signal_power_watts,noise_power_watts,snr_db] = ddm_snr2watts(ddm_snr_db,noise_std_ch3,cal_factors_ch3);

end
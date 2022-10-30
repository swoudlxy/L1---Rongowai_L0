% this function uses "old" CYGNSS way to derive noise floor from individual
% DDMs

function [signal_power_watts,noise_power_watts,noise_floor_counts] = L1a_counts2watts1(raw_counts,ANZ_port)

noise_floor_counts = mean(raw_counts(:,end-4:end),'all');
ddm_signal_counts = raw_counts-noise_floor_counts;

signal_cal_factors_db = [-162.93,-163.98,-161.24];      % calibration factors
cable_loss_db = [1.8051,0.6600,0.5840];                 % cable loss

% perform calibration based on the input ANZ port channel
switch ANZ_port

    case 1
        signal_cal_factors_db_ch = signal_cal_factors_db(1);
        cable_loss_db_ch = cable_loss_db(1);
        
    case 2
        signal_cal_factors_db_ch = signal_cal_factors_db(2);
        cable_loss_db_ch = cable_loss_db(2);        
            
    case 3
        signal_cal_factors_db_ch = signal_cal_factors_db(3);
        cable_loss_db_ch = cable_loss_db(3);

end

signal_power_watts = ddm_signal_counts*db2pow(signal_cal_factors_db_ch)*db2pow(cable_loss_db_ch);
noise_power_watts = noise_floor_counts*db2pow(signal_cal_factors_db_ch)*db2pow(cable_loss_db_ch);

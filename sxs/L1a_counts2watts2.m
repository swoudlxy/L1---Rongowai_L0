% this function calibrates from raw DDM counts to DDM power in watts

function ddm_power_watts = L1a_counts2watts2(ddm_counts,ANZ_port, ...
    ddm_counts_cal,ddm_power_cal,std_dev,noise_floor)

binning_thres_db = [50.5 49.6 50.4];
binning_thres = db2pow(binning_thres_db);

noise_floor_LHCP = noise_floor(1);
noise_floor_RHCP = noise_floor(2);

% select approiate calibration constants based on the input ANZ port channel
switch ANZ_port

    case 1      % no case 1 in real case because this is zenith channel
        ddm_counts_ch = ddm_counts_cal(1,:);
        ddm_power_ch = ddm_power_cal(1,:);
        
        mag_std_dev_ch = std_dev(1);

        binning_thres_ch = binning_thres(1);
        
    case 2
        ddm_counts_ch = ddm_counts_cal(2,:);
        ddm_power_ch = ddm_power_cal(2,:);

        ddm_counts = ddm_counts-noise_floor_LHCP;

        mag_std_dev_ch = std_dev(2);

        binning_thres_ch = binning_thres(2); 
            
    case 3
        ddm_counts_ch = ddm_counts_cal(3,:);
        ddm_power_ch = ddm_power_cal(3,:);

        ddm_counts = ddm_counts-noise_floor_RHCP;

        mag_std_dev_ch = std_dev(3);

        binning_thres_ch = binning_thres(3);

end

std_dev_ch = mag_std_dev_ch^2;              % standard deviation also needs to be in power scale
ddm_power = zeros(size(ddm_counts));

% evaluate ddm power in watts
ddm_power(:) = interp1(ddm_counts_ch,ddm_power_ch,ddm_counts(:),'spline');
ddm_power_watts = ddm_power*std_dev_ch/binning_thres_ch;                  
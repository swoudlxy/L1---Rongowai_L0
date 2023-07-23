% this function calibrates from raw DDM counts to DDM power in watts

function ddm_power_watts = L1a_counts2watts(ddm_counts,ANZ_port,ddm_counts_cal_db,ddm_power_cal_dbm,std_dev)

binning_thres_db = [50.5 49.6 50.4];

% select approiate calibration constants based on the input ANZ port channel
switch ANZ_port

    case 1
        ddm_counts_db_ch = ddm_counts_cal_db(1,:);
        ddm_power_dbm_ch = ddm_power_cal_dbm(1,:);
        
        std_dev_ch = std_dev(1);

        binning_thres_db_ch = binning_thres_db(1);
        
    case 2
        ddm_counts_db_ch = ddm_counts_cal_db(2,:);
        ddm_power_dbm_ch = ddm_power_cal_dbm(2,:)-11.1;     % LHCP minus 11.1 dB - 27 June

        std_dev_ch = std_dev(2);

        binning_thres_db_ch = binning_thres_db(2); 
            
    case 3
        ddm_counts_db_ch = ddm_counts_cal_db(3,:);
        ddm_power_dbm_ch = ddm_power_cal_dbm(3,:)-14.4;     % RHCP minus 14.4 dB - 27 June

        std_dev_ch = std_dev(3);

        binning_thres_db_ch = binning_thres_db(3);

end

% convert to dB scale
ddm_counts_db = pow2db(ddm_counts);
std_dev_db_ch = mag2db(std_dev_ch);

ddm_power_dbm = zeros(size(ddm_counts_db));

% evaluate ddm power in dBm
ddm_power_dbm(:) = interp1(ddm_counts_db_ch,ddm_power_dbm_ch,ddm_counts_db(:),'spline');
ddm_power_dbm = ddm_power_dbm+std_dev_db_ch-binning_thres_db_ch;    % cable loss to be compensated when computing BRCS

% convert to watts
% remove 9 dB difference which has been compensated above 27 June
ddm_power_watts = db2pow(ddm_power_dbm-30);                         
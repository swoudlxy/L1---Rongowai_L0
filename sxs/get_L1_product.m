% This function solves all L1 variables and packets as a structure for
% processing multiple L0 files

% L1 Algorithm version: 1.2.1
% L1 data version: 1.2.0

function L1_postCal = get_L1_product(   L0_filename, ...
                                        L1a_cal_ddm_counts_db,L1a_cal_ddm_power_dbm, ...
                                        dem,dtu10,landmask_nz,lcv_mask,water_mask, ...
                                        SV_PRN_LUT,SV_eirp_LUT,LHCP_pattern,RHCP_pattern, ...
                                        phy_ele_size)

% Prelaunch 1: Load L0 data and filter invalid values

% PVT GPS week and sec
pvt_gps_week = double(ncread(L0_filename,'/science/GPS_week_of_SC_attitude'));
pvt_gps_sec = double(ncread(L0_filename,'/science/GPS_second_of_SC_attitude'));

% rx positions in ECEF, metres
rx_pos_x_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_position_x_ecef_m'));
rx_pos_y_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_position_y_ecef_m'));
rx_pos_z_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_position_z_ecef_m'));

% rx velocity in ECEF, m/s
rx_vel_x_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_velocity_x_ecef_mps'));
rx_vel_y_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_velocity_y_ecef_mps'));
rx_vel_z_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_velocity_z_ecef_mps'));

% rx attitude, deg
rx_pitch_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_attitude_pitch_deg'));
rx_roll_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_attitude_roll_deg'));
rx_yaw_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_attitude_yaw_deg'));

% rx clock bias and drifts
rx_clk_bias_m_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_clock_bias_m'));
rx_clk_drift_mps_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_clock_drift_mps'));

% tx ID/satellite PRN
transmitter_id = double(ncread(L0_filename,'/science/ddm/transmitter_id'));

% raw counts and ddm parameters
first_scale_factor = double(ncread(L0_filename,'/science/ddm/first_scale_factor'));
raw_counts = double(ncread(L0_filename,'/science/ddm/counts'));                         % raw counts, uncalibrated
zenith_i2q2 = double(ncread(L0_filename,'/science/ddm/zenith_i2_plus_q2'));             % zenith counts

rf_source = double(ncread(L0_filename,'/science/ddm/RF_source'));                       % RF source

% binning standard deviation
std_dev_rf1 = double(ncread(L0_filename,'/science/ddm/RF1_zenith_RHCP_std_dev'));
std_dev_rf2 = double(ncread(L0_filename,'/science/ddm/RF2_nadir_LHCP_std_dev')); 
std_dev_rf3 = double(ncread(L0_filename,'/science/ddm/RF3_nadir_RHCP_std_dev'));

% absolute ddm center delay and doppler
delay_center_chips = double(ncread(L0_filename,'/science/ddm/center_delay_bin_code_phase'));    
doppler_center_hz = double(ncread(L0_filename,'/science/ddm/center_doppler_bin_frequency'));

% coherent duration and noncoherent integration
coherent_duration = double(ncread(L0_filename,'/science/ddm/L1_E1_coherent_duration'));
non_coherent_integrations = double(ncread(L0_filename,'/science/ddm/L1_E1_non_coherent_integrations'));

% NGRx estimate additional delay path
add_range_to_sp_pvt = double(ncread(L0_filename,'/science/ddm/additional_range_to_SP'));

% antenna temperatures and engineering timestamp
eng_timestamp = double(ncread(L0_filename,'/eng/packet_creation_time'));
zenith_ant_temp_eng = double(ncread(L0_filename,'/eng/zenith_ant_temp'));
nadir_ant_temp_eng = double(ncread(L0_filename,'/eng/nadir_ant_temp'));

% the below scrpits filter valid data

% rx-related variables
index1 = pvt_gps_week>0;

pvt_gps_week = pvt_gps_week(index1);        pvt_gps_sec = pvt_gps_sec(index1);

rx_pos_x_pvt = rx_pos_x_pvt(index1);        rx_pos_y_pvt = rx_pos_y_pvt(index1);
rx_pos_z_pvt = rx_pos_z_pvt(index1);

rx_vel_x_pvt = rx_vel_x_pvt(index1);        rx_vel_y_pvt = rx_vel_y_pvt(index1);
rx_vel_z_pvt = rx_vel_z_pvt(index1);

rx_roll_pvt = rx_roll_pvt(index1);          rx_pitch_pvt = rx_pitch_pvt(index1);
rx_yaw_pvt = rx_yaw_pvt(index1);

rx_clk_bias_m_pvt = rx_clk_bias_m_pvt(index1);
rx_clk_drift_mps_pvt = rx_clk_drift_mps_pvt(index1);

% remove first and last few zeros
index0_b = find(pvt_gps_week>0,1);
index0_e = find(pvt_gps_week>0,1,'last');

pvt_gps_week = pvt_gps_week(index0_b:index0_e);
pvt_gps_sec = pvt_gps_sec(index0_b:index0_e);

rx_pos_x_pvt = rx_pos_x_pvt(index0_b:index0_e);
rx_pos_y_pvt = rx_pos_y_pvt(index0_b:index0_e);
rx_pos_z_pvt = rx_pos_z_pvt(index0_b:index0_e);

rx_vel_x_pvt = rx_vel_x_pvt(index0_b:index0_e);
rx_vel_y_pvt = rx_vel_y_pvt(index0_b:index0_e);
rx_vel_z_pvt = rx_vel_z_pvt(index0_b:index0_e);

rx_roll_pvt = rx_roll_pvt(index0_b:index0_e);
rx_pitch_pvt = rx_pitch_pvt(index0_b:index0_e);
rx_yaw_pvt = rx_yaw_pvt(index0_b:index0_e);

rx_clk_bias_m_pvt = rx_clk_bias_m_pvt(index0_b:index0_e);
rx_clk_drift_mps_pvt = rx_clk_drift_mps_pvt(index0_b:index0_e);

% identify and compensate the value equal to 0 (randomly happens)
index0 = find(pvt_gps_week==0);
if ~isempty(index0)
    L = length(index0);

    for l = 1:L

        index0_1 = index0(l);

        pvt_gps_week(index0_1) = mean([pvt_gps_week(index0_1-1) pvt_gps_week(index0_1+1)]);
        pvt_gps_sec(index0_1) = mean([pvt_gps_sec(index0_1-1) pvt_gps_sec(index0_1+1)]);

        rx_pos_x_pvt(index0_1) = mean([rx_pos_x_pvt(index0_1-1) rx_pos_x_pvt(index0_1+1)]);
        rx_pos_y_pvt(index0_1) = mean([rx_pos_y_pvt(index0_1-1) rx_pos_y_pvt(index0_1+1)]);
        rx_pos_z_pvt(index0_1) = mean([rx_pos_z_pvt(index0_1-1) rx_pos_z_pvt(index0_1+1)]);

        rx_vel_x_pvt(index0_1) = mean([rx_vel_x_pvt(index0_1-1) rx_vel_x_pvt(index0_1+1)]);
        rx_vel_y_pvt(index0_1) = mean([rx_vel_y_pvt(index0_1-1) rx_vel_y_pvt(index0_1+1)]);
        rx_vel_z_pvt(index0_1) = mean([rx_vel_z_pvt(index0_1-1) rx_vel_z_pvt(index0_1+1)]);

        rx_roll_pvt(index0_1) = mean([rx_roll_pvt(index0_1-1) rx_roll_pvt(index0_1+1)]);
        rx_pitch_pvt(index0_1) = mean([rx_pitch_pvt(index0_1-1) rx_pitch_pvt(index0_1+1)]);
        rx_yaw_pvt(index0_1) = mean([rx_yaw_pvt(index0_1-1) rx_yaw_pvt(index0_1+1)]);

        rx_clk_bias_m_pvt(index0_1) = mean([rx_clk_bias_m_pvt(index0_1-1) rx_clk_bias_m_pvt(index0_1+1)]);
        rx_clk_drift_mps_pvt(index0_1) = mean([rx_clk_drift_mps_pvt(index0_1-1) rx_clk_drift_mps_pvt(index0_1+1)]);

    end

end

% ddm-related variables
index2 = ~isnan(transmitter_id(1,:));

transmitter_id = transmitter_id(:,index2);

first_scale_factor = first_scale_factor(:,index2);
raw_counts = raw_counts(:,:,:,index2); 
zenith_i2q2 = zenith_i2q2(:,index2);

rf_source = rf_source(:,index2);

std_dev_rf1 = std_dev_rf1(index2);
std_dev_rf2 = std_dev_rf2(index2);
std_dev_rf3 = std_dev_rf3(index2);

% absolute ddm center delay and doppler
delay_center_chips = delay_center_chips(:,index2);
doppler_center_hz = doppler_center_hz(:,index2);

% coherent duration and noncoherent integration
coherent_duration = coherent_duration(index2)/1000;                     % convert to seconds
non_coherent_integrations = non_coherent_integrations(index2)/1000;

% NGRx estimate additional delay path
add_range_to_sp_pvt = add_range_to_sp_pvt(:,index2);

% the below is to process when ddm-related and rx-related variables do not
% have the same length, which happens for some of the L0 products
diff = length(pvt_gps_week)-length(transmitter_id);

if diff > 0

    offset_idx = diff+1;

    pvt_gps_week = pvt_gps_week(offset_idx:end);
    pvt_gps_sec = pvt_gps_sec(offset_idx:end);

    rx_pos_x_pvt = rx_pos_x_pvt(offset_idx:end); 
    rx_pos_y_pvt = rx_pos_y_pvt(offset_idx:end);
    rx_pos_z_pvt = rx_pos_z_pvt(offset_idx:end);

    rx_vel_x_pvt = rx_vel_x_pvt(offset_idx:end); 
    rx_vel_y_pvt = rx_vel_y_pvt(offset_idx:end);
    rx_vel_z_pvt = rx_vel_z_pvt(offset_idx:end);

    rx_roll_pvt = rx_roll_pvt(offset_idx:end);   
    rx_pitch_pvt = rx_pitch_pvt(offset_idx:end);
    rx_yaw_pvt = rx_yaw_pvt(offset_idx:end);

    rx_clk_bias_m_pvt = rx_clk_bias_m_pvt(offset_idx:end);
    rx_clk_drift_mps_pvt = rx_clk_drift_mps_pvt(offset_idx:end);

elseif diff < 0

    offset_idx = abs(diff)+1;

    transmitter_id = transmitter_id(:,offset_idx:end);

    first_scale_factor = first_scale_factor(:,offset_idx:end);
    raw_counts = raw_counts(:,:,:,offset_idx:end); 
    zenith_i2q2 = zenith_i2q2(:,offset_idx:end);

    rf_source = rf_source(:,offset_idx:end);

    std_dev_rf1 = std_dev_rf1(offset_idx:end);
    std_dev_rf2 = std_dev_rf2(offset_idx:end);
    std_dev_rf3 = std_dev_rf3(offset_idx:end);
        
    delay_center_chips = delay_center_chips(:,offset_idx:end);
    doppler_center_hz = doppler_center_hz(:,offset_idx:end);

    % coherent duration and noncoherent integration
    coherent_duration = coherent_duration(offset_idx:end)/1000;
    non_coherent_integrations = non_coherent_integrations(offset_idx:end)/1000;

    % NGRx estimate additional delay path
    add_range_to_sp_pvt = add_range_to_sp_pvt(:,offset_idx:end);

end

% temperatures from engineering data
index3 = ~isnan(eng_timestamp);

eng_timestamp = eng_timestamp(index3);
nadir_ant_temp_eng = nadir_ant_temp_eng(index3);
zenith_ant_temp_eng = zenith_ant_temp_eng(index3);

invalid = nan;                              % defines the value to be used for invalid fields
I = length(pvt_gps_sec);                    % total length of samples
J = 20;                                     % maximal NGRx capacity

% initialise a structure to save L1 results
L1_postCal = struct;

% Part 1: General processing
% This part derives global constants, timestamps, and all the other
% parameters at ddm timestamps

% initialise output data array for Part 1 processing
pvt_utc = zeros(I,1)+invalid;       ddm_utc = zeros(I,1)+invalid;
gps_week = zeros(I,1)+invalid;      gps_tow = zeros(I,1)+invalid;
ddm_pvt_bias = zeros(I,1)+invalid;
add_range_to_sp = zeros(J,I)+invalid;
status_flags_one_hz = zeros(I,1)+invalid;

% derive and save ddm_timestamp_utc/gps to L1 structure
for i = 1:I

    pvt_gps_week1 = pvt_gps_week(i);
    pvt_gps_sec1 = pvt_gps_sec(i);
    non_coherent1 = non_coherent_integrations(i);

    % convert pvt_gps_time to pvt_utc time
    D_pvt1 = gpstime2utc(pvt_gps_week1,pvt_gps_sec1);
    D_pvt2 = datetime(D_pvt1(1),D_pvt1(2),D_pvt1(3),D_pvt1(4),D_pvt1(5),D_pvt1(6)); % human time
    pvt_utc1 = convertTo(D_pvt2,'posixtime');    

    % derive ddm_utc, mid point of the non-coherent integrations
    ddm_utc1 = pvt_utc1+non_coherent1/2;    
    
    % derive ddm_gps_week and ddm_gps_sec
    D_ddm1 = datetime(ddm_utc1,'ConvertFrom','posixtime');
    [year,month,date] = ymd(D_ddm1);
    [hour,mins,secs] = hms(D_ddm1);
    [ddm_gps_week1,ddm_gps_sec1] = utc2gpstime(year,month,date,hour,mins,secs);

    % save pvt and ddm timestamps for interpolation
    pvt_utc(i) = pvt_utc1;          ddm_utc(i) = ddm_utc1;
    gps_week(i) = ddm_gps_week1;    gps_tow(i) = ddm_gps_sec1;
    ddm_pvt_bias(i) = non_coherent1/2;

end

% linear interpolation all the values at ddm timestamp
rx_pos_x = interp1(pvt_utc,rx_pos_x_pvt,ddm_utc,'linear','extrap');
rx_pos_y = interp1(pvt_utc,rx_pos_y_pvt,ddm_utc,'linear','extrap');
rx_pos_z = interp1(pvt_utc,rx_pos_z_pvt,ddm_utc,'linear','extrap');
rx_pos_xyz = [rx_pos_x rx_pos_y rx_pos_z];

rx_vel_x = interp1(pvt_utc,rx_vel_x_pvt,ddm_utc,'linear','extrap');
rx_vel_y = interp1(pvt_utc,rx_vel_y_pvt,ddm_utc,'linear','extrap');
rx_vel_z = interp1(pvt_utc,rx_vel_z_pvt,ddm_utc,'linear','extrap');
rx_vel_xyz = [rx_vel_x rx_vel_y rx_vel_z];

rx_roll = interp1(pvt_utc,rx_roll_pvt,ddm_utc,'linear','extrap');
rx_pitch = interp1(pvt_utc,rx_pitch_pvt,ddm_utc,'linear','extrap');
rx_yaw = interp1(pvt_utc,rx_yaw_pvt,ddm_utc,'linear','extrap');
rx_attitude = [rx_roll rx_pitch rx_yaw];

rx_clk_bias_m = interp1(pvt_utc,rx_clk_bias_m_pvt,ddm_utc,'linear','extrap');
rx_clk_drift_mps = interp1(pvt_utc,rx_clk_drift_mps_pvt,ddm_utc,'linear','extrap');

for j = 1:J
    
    temp = add_range_to_sp_pvt(j,:);
    temp1 = interp1(pvt_utc,temp,ddm_utc,'linear','extrap');
    add_range_to_sp(j,:) = temp1;

end

ant_temp_zenith = interp1(eng_timestamp,zenith_ant_temp_eng,ddm_utc,'linear','extrap');
ant_temp_nadir = interp1(eng_timestamp,nadir_ant_temp_eng,ddm_utc,'linear','extrap');

% get <lat,lon,alt> of aircraft and status flag
rx_pos_lla = ecef2lla(rx_pos_xyz);
for i = 1:I

    dist1 = get_map_value(rx_pos_lla(i,1),rx_pos_lla(i,2),landmask_nz);

    if dist1 <= 0
        status_flags_one_hz(i) = 4;
    elseif dist1 > 0
        status_flags_one_hz(i) = 5;
    end

end

% write global variables
time_coverage_start = string(datetime(ddm_utc(1),'ConvertFrom','posixtime'));
L1_postCal.time_coverage_start = time_coverage_start;

time_coverage_end = string(datetime(ddm_utc(end),'ConvertFrom','posixtime'));
L1_postCal.time_coverage_end = time_coverage_end;

L1_postCal.time_coverage_resolution = ddm_utc(2)-ddm_utc(1);

% time coverage
time_duration = ddm_utc(end)-ddm_utc(1)+1;
hours = floor(time_duration/3600);
minutes = floor((time_duration-hours*3600)/60);
seconds = time_duration-hours*3600-minutes*60;
time_coverage_duration = ['P0DT' num2str(hours) 'H' num2str(minutes) 'M' num2str(seconds) 'S'];

L1_postCal.time_coverage_duration = time_coverage_duration;

% below is new for algorithm version 1.1
ref_timestamp_utc = ddm_utc(1);

pvt_timestamp_utc = pvt_utc-ref_timestamp_utc;
ddm_timestamp_utc = ddm_utc-ref_timestamp_utc;
% ends

L1_postCal.aircraft_reg = 'ZK-NFA';             % default value
L1_postCal.ddm_source = 2;                      % 1 = GPS signal simulator, 2 = aircraft
L1_postCal.ddm_time_type_selector = 1;          % 1 = middle of DDM sampling period
L1_postCal.delay_resolution = 0.25;             % unit in chips
L1_postCal.dopp_resolution = 500;               % unit in Hz
L1_postCal.dem_source = 'SRTM30';

% write algorithm and LUT versions
L1_postCal.l1_algorithm_version = '1.1';
L1_postCal.l1_data_version = '1.12';
L1_postCal.l1a_sig_LUT_version = '1';
L1_postCal.l1a_noise_LUT_version = '1';
L1_postCal.ngrx_port_mapping_version = '1';
L1_postCal.nadir_ant_data_version = '1';
L1_postCal.zenith_ant_data_version = '1';
L1_postCal.prn_sv_maps_version = '1';
L1_postCal.gps_eirp_param_version = '7';
L1_postCal.land_mask_version = '1';
L1_postCal.surface_type_version = '1';
L1_postCal.mean_sea_surface_version = '1';
L1_postCal.per_bin_ant_version = '1';

% write timestamps and ac-related variables
L1_postCal.pvt_timestamp_gps_week = pvt_gps_week;
L1_postCal.pvt_timestamp_gps_sec = pvt_gps_sec;
L1_postCal.pvt_timestamp_utc = pvt_timestamp_utc;       % algorithm version 1.1

L1_postCal.ddm_timestamp_gps_week = gps_week;
L1_postCal.ddm_timestamp_gps_sec = gps_tow;
L1_postCal.ddm_timestamp_utc = ddm_timestamp_utc;       % algorithm version 1.1

L1_postCal.ddm_pvt_bias = ddm_pvt_bias;

% 0-indexed sample and DDM
L1_postCal.sample = (0:1:I-1)';
L1_postCal.ddm = (0:1:J-1)';

L1_postCal.sp_fsw_delay = delay_center_chips;
L1_postCal.sp_ngrx_dopp = doppler_center_hz;

L1_postCal.add_range_to_sp = add_range_to_sp;
L1_postCal.add_range_to_sp_pvt = add_range_to_sp_pvt;

L1_postCal.ac_lat = rx_pos_lla(:,1);
L1_postCal.ac_lon = rx_pos_lla(:,2);
L1_postCal.ac_alt = rx_pos_lla(:,3);

L1_postCal.ac_pos_x_pvt = rx_pos_x_pvt;
L1_postCal.ac_pos_y_pvt = rx_pos_y_pvt;
L1_postCal.ac_pos_z_pvt = rx_pos_z_pvt;

L1_postCal.ac_pos_x = rx_pos_x;
L1_postCal.ac_pos_y = rx_pos_y;
L1_postCal.ac_pos_z = rx_pos_z;

L1_postCal.ac_vel_x_pvt = rx_vel_x_pvt;
L1_postCal.ac_vel_y_pvt = rx_vel_y_pvt;
L1_postCal.ac_vel_z_pvt = rx_vel_z_pvt;

L1_postCal.ac_vel_x = rx_vel_x;
L1_postCal.ac_vel_y = rx_vel_y;
L1_postCal.ac_vel_z = rx_vel_z;

L1_postCal.ac_roll_pvt = rx_roll_pvt;
L1_postCal.ac_pitch_pvt = rx_pitch_pvt;
L1_postCal.ac_yaw_pvt = rx_yaw_pvt;

L1_postCal.ac_roll = rx_attitude(:,1);
L1_postCal.ac_pitch = rx_attitude(:,2);
L1_postCal.ac_yaw = rx_attitude(:,3);

L1_postCal.rx_clk_bias_pvt = rx_clk_bias_m_pvt;
L1_postCal.rx_clk_drift_pvt = rx_clk_drift_mps_pvt;

L1_postCal.rx_clk_bias = rx_clk_bias_m;
L1_postCal.rx_clk_drift = rx_clk_drift_mps;

L1_postCal.ant_temp_nadir = ant_temp_nadir;
L1_postCal.ant_temp_zenith = ant_temp_zenith;

L1_postCal.status_flags_one_hz = status_flags_one_hz;

% Part 1 ends

% Part 2: Derive TX related variables
% This part derives TX positions and velocities, maps between PRN and SVN,
% and gets track ID

trans_id_unique = unique(transmitter_id);
trans_id_unique = trans_id_unique(trans_id_unique>0);

% new SP3 naming policy
D_ddm1 = datetime(ddm_utc(1),'ConvertFrom','posixtime');
[year1,~,~] = ymd(D_ddm1);
doy1 = day(D_ddm1,'dayofyear');

D_ddm2 = datetime(ddm_utc(end),'ConvertFrom','posixtime');
[year2,~,~] = ymd(D_ddm1);
doy2 = day(D_ddm2,'dayofyear');

if doy1 == doy2
    flag = 0;
else
    flag = 1;
    dow = floor(gps_tow/86400);         % day of week
    change_idx = find(ischange(dow)==1);
end

% compensate 0s in front of doy
if doy1 < 10
    doy1 = ['00' num2str(doy1)];
elseif doy1 < 100
    doy1 = ['0' num2str(doy1)];
else
    doy1 = num2str(doy1);
end

if doy2 < 10
    doy2 = ['00' num2str(doy2)];
elseif doy2 < 100
    doy2 = ['0' num2str(doy2)];
else
    doy2 = num2str(doy2);
end

% initalise variables 
tx_pos_x = zeros(J,I)+invalid;      tx_vel_x = zeros(J,I)+invalid;
tx_pos_y = zeros(J,I)+invalid;      tx_vel_y = zeros(J,I)+invalid;
tx_pos_z = zeros(J,I)+invalid;      tx_vel_z = zeros(J,I)+invalid;

tx_clk_bias = zeros(J,I)+invalid;

prn_code = zeros(J,I)+invalid;      sv_num = zeros(J,I)+invalid;
track_id = zeros(J,I)+invalid;

for i = 1:I

    % gps timestamp
    ddm_gps_timestamp1.gps_week = gps_week(i); 
    ddm_gps_timestamp1.gps_tow = gps_tow(i);

    for j = 1:J/2
    
        % assign PRN
        transmitter_id1 = transmitter_id(j,i);
        prn1 = transmitter_id1;    

        if prn1 ~= 0

            sv_num1 = SV_PRN_LUT(SV_PRN_LUT(:,1)==prn1,2);  % assign SVN

            % assign correct sp3 files if the flight cross two dates
            if flag == 0

                sp3_filename = ['IGS0OPSRAP_' num2str(year1) doy1 '0000_01D_15M_ORB.SP3'];
                gps_orbit_filename = ['..//dat//orbits//' sp3_filename];                

            elseif flag == 1
                if i<change_idx
                    sp3_filename = ['IGS0OPSRAP_' num2str(year1) doy1 '0000_01D_15M_ORB.SP3'];
                    gps_orbit_filename = ['..//dat//orbits//' sp3_filename];

                elseif i>=change_idx
                    sp3_filename = ['IGS0OPSRAP_' num2str(year2) doy2 '0000_01D_15M_ORB.SP3'];
                    gps_orbit_filename = ['..//dat//orbits//' sp3_filename];
                
                end

            end

            [tx_pos_xyz1,tx_vel_xyz1,tx_clk_bias1,~] = gps_posvel(prn1,ddm_gps_timestamp1, ...
                gps_orbit_filename);

            tx_pos_x(j,i) = tx_pos_xyz1(1); tx_vel_x(j,i) = tx_vel_xyz1(1);
            tx_pos_y(j,i) = tx_pos_xyz1(2); tx_vel_y(j,i) = tx_vel_xyz1(2);
            tx_pos_z(j,i) = tx_pos_xyz1(3); tx_vel_z(j,i) = tx_vel_xyz1(3);

            tx_clk_bias(j,i) = tx_clk_bias1;           

            prn_code(j,i) = prn1;           
            sv_num(j,i) = sv_num1;
            track_id(j,i) = find(trans_id_unique == transmitter_id1);

        end

    end
end

%{
% old sp3 naming policy
gps_week1 = num2str(gps_week(1));
gps_dow1 = num2str(floor(gps_tow(1)/86400));

gps_week2 = num2str(gps_week(end));
gps_dow2 = num2str(floor(gps_tow(end)/86400));

if dow1 == dow2
    flag = 0;
else
    flag = 1;
    dow = floor(gps_tow/86400);         % day of week
    change_idx = find(ischange(dow)==1);
end

% initalise variables 
tx_pos_x = zeros(J,I)+invalid;      tx_vel_x = zeros(J,I)+invalid;
tx_pos_y = zeros(J,I)+invalid;      tx_vel_y = zeros(J,I)+invalid;
tx_pos_z = zeros(J,I)+invalid;      tx_vel_z = zeros(J,I)+invalid;

tx_clk_bias = zeros(J,I)+invalid;

prn_code = zeros(J,I)+invalid;      sv_num = zeros(J,I)+invalid;
track_id = zeros(J,I)+invalid;

for i = 1:I

    % gps timestamp
    ddm_gps_timestamp1.gps_week = gps_week(i); 
    ddm_gps_timestamp1.gps_tow = gps_tow(i);

    for j = 1:J/2

        % assign PRN
        transmitter_id1 = transmitter_id(j,i);
        prn1 = transmitter_id1;    

        if prn1 ~= 0

            sv_num1 = SV_PRN_LUT(SV_PRN_LUT(:,1)==prn1,2);  % assign SVN

            % assign correct sp3 files if the flight cross two dates
            if flag == 0

                sp3_filename = ['igr' num2str(gps_week1) num2str(gps_dow1) '.sp3'];
                gps_orbit_filename = ['..//dat//orbits//' sp3_filename];                

            elseif flag == 1
                if i<change_idx
                    sp3_filename = ['igr' num2str(gps_week1) num2str(gps_dow1) '.sp3'];
                    gps_orbit_filename = ['..//dat//orbits//' sp3_filename];

                elseif i>=change_idx
                    sp3_filename = ['igr' num2str(gps_week2) num2str(gps_dow2) '.sp3'];
                    gps_orbit_filename = ['..//dat//orbits//' sp3_filename];
                
                end

            end

            [tx_pos_xyz1,tx_vel_xyz1,tx_clk_bias1,~] = gps_posvel(prn1,ddm_gps_timestamp1, ...
                gps_orbit_filename);

            tx_pos_x(j,i) = tx_pos_xyz1(1); tx_vel_x(j,i) = tx_vel_xyz1(1);
            tx_pos_y(j,i) = tx_pos_xyz1(2); tx_vel_y(j,i) = tx_vel_xyz1(2);
            tx_pos_z(j,i) = tx_pos_xyz1(3); tx_vel_z(j,i) = tx_vel_xyz1(3);

            tx_clk_bias(j,i) = tx_clk_bias1;
           

            prn_code(j,i) = prn1;           
            sv_num(j,i) = sv_num1;
            track_id(j,i) = find(trans_id_unique == transmitter_id1);
    
        end

    end
end
%}

% extend to RHCP channels
tx_pos_x(J/2+1:J,:) = tx_pos_x(1:J/2,:);    tx_vel_x(J/2+1:J,:) = tx_vel_x(1:J/2,:);
tx_pos_y(J/2+1:J,:) = tx_pos_y(1:J/2,:);    tx_vel_x(J/2+1:J,:) = tx_vel_x(1:J/2,:);
tx_pos_z(J/2+1:J,:) = tx_pos_z(1:J/2,:);    tx_vel_x(J/2+1:J,:) = tx_vel_x(1:J/2,:);

tx_clk_bias(J/2+1:J,:) = tx_clk_bias(1:J/2,:);

prn_code(J/2+1:J,:) = prn_code(1:J/2,:);
sv_num(J/2+1:J,:) = sv_num(1:J/2,:);
track_id(J/2+1:J,:) = track_id(1:J/2,:);

% write TX variables
L1_postCal.tx_pos_x = tx_pos_x;
L1_postCal.tx_pos_y = tx_pos_y;
L1_postCal.tx_pos_z = tx_pos_z;

L1_postCal.tx_vel_x = tx_vel_x;
L1_postCal.tx_vel_y = tx_vel_y;
L1_postCal.tx_vel_z = tx_vel_z;

L1_postCal.tx_clk_bias = tx_clk_bias;

L1_postCal.prn_code = prn_code;
L1_postCal.sv_num = sv_num;
L1_postCal.track_id = track_id;

% Part 3: L1a calibration
% this part converts from raw counts to signal power in watts and complete
% L1a calibration

offset = 4;                         % offset delay rows to derive noise floor

% initialise variables for L1a results
ddm_power_counts = zeros(5,40,J,I)+invalid;
power_analog = zeros(5,40,J,I)+invalid;

%noise_floor_counts = zeros(J,I)+invalid;
%noise_floor = zeros(J,I)+invalid;
snr_db = zeros(J,I)+invalid;

peak_ddm_counts = zeros(J,I)+invalid;
peak_ddm_watts = zeros(J,I)+invalid;
peak_delay_bin = zeros(J,I)+invalid;

ddm_noise_counts = zeros(J,I)+invalid;
%ddm_noise_watts = zeros(J,I)+invalid;

ddm_ant = zeros(J,I)+invalid;
%inst_gain = zeros(J,I)+invalid;

% derive signal power
for i = 1:I

    % retrieve noise standard deviation in counts for all three channels
    std_dev1 = [std_dev_rf1(i),std_dev_rf2(i),std_dev_rf3(i)];
    
    for j = 1:J

        prn_code1 = prn_code(j,i);
        rf_source1 = rf_source(j,i);

        first_scale_factor1 = first_scale_factor(j,i);
        raw_counts1 = raw_counts(:,:,j,i);

        % solve only when presenting a valid PRN and DDM counts
        if      (~isnan(prn_code1)) && ...
                (raw_counts1(1,1) ~= raw_counts1(3,21)) && ...
                (raw_counts1(1,1) ~=0)

            ANZ_port1 = get_ANZ_port(rf_source1);
            ddm_power_counts1 = raw_counts1*first_scale_factor1;

            % perform L1a calibration from Counts to Watts
            ddm_power_watts1 = L1a_counts2watts(ddm_power_counts1,ANZ_port1, ...
                L1a_cal_ddm_counts_db,L1a_cal_ddm_power_dbm,std_dev1);

            % noise floor in counts for each DDM
            % updated for ver. 1.21
            ddm_noise_counts1 = mean(ddm_power_counts1(:,end-offset:end),'all');
            %ddm_noise_watts1 = L1a_counts2watts(ddm_noise_counts1,ANZ_port1, ...
            %    L1a_cal_ddm_counts_db,L1a_cal_ddm_power_dbm,std_dev1);

            % peak ddm location
            peak_ddm_counts1 = max(max(ddm_power_counts1));
            [~,peak_delay_bin1] = find(ddm_power_counts1==peak_ddm_counts1,1);

            peak_ddm_watts1 = max(max(ddm_power_watts1));
            
            % save variables
            ddm_power_counts(:,:,j,i) = ddm_power_counts1;
            power_analog(:,:,j,i) = ddm_power_watts1;
            
            ddm_ant(j,i) = ANZ_port1;            

            ddm_noise_counts(j,i) = ddm_noise_counts1;
            %ddm_noise_watts(j,i) = ddm_noise_watts1;

            peak_ddm_counts(j,i) = peak_ddm_counts1;
            peak_ddm_watts(j,i) = peak_ddm_watts1;
            peak_delay_bin(j,i) = peak_delay_bin1;

        end
        
    end    
end

% changed in ver 1.21, may bring back in a future version
% noise floor, SNR and instrument gain will be solved after getting
% specular coordinates
%{
% derive noise floor, SNR and instrument gain
for i = 1:I

    noise_counts_LHCP1 = ddm_noise_counts(1:J/2,i);
    noise_watts_LHCP1 = ddm_noise_watts(1:J/2,i);
    peak_delay_bin_LHCP1 = peak_delay_bin(1:J/2,i);

    noise_counts_RHCP1 = ddm_noise_counts(J/2+1:J,i);
    noise_watts_RHCP1 = ddm_noise_watts(J/2+1:J,i);
    peak_delay_bin_RHCP1 = peak_delay_bin(J/2+1:J,i);

    noise_index_LHCP = find(peak_delay_bin_LHCP1<31 & peak_delay_bin_LHCP1>0);
    noise_index_RHCP = find(peak_delay_bin_RHCP1<31 & peak_delay_bin_RHCP1>0);

    if ~isempty(noise_index_LHCP)
        avg_noise_counts_LHCP1 = mean(noise_counts_LHCP1(noise_index_LHCP));
        avg_noise_watts_LHCP1 = mean(noise_watts_LHCP1(noise_index_LHCP));

        avg_noise_counts_RHCP1 = mean(noise_counts_RHCP1(noise_index_RHCP));
        avg_noise_watts_RHCP1 = mean(noise_watts_RHCP1(noise_index_RHCP));

    elseif (isempty(noise_index_LHCP)) && (sum(~isnan(noise_counts_LHCP1))>0)
        avg_noise_counts_LHCP1 = noise_floor_counts(j,i-1);
        avg_noise_watts_LHCP1 = noise_floor(j,i-1);

        avg_noise_counts_RHCP1 = noise_floor_counts(j+J/2,i-1);
        avg_noise_watts_RHCP1 = noise_floor(j+J/2,i-1);

    end

    for j = 1:J/2

        peak_power_counts_LHCP1 = peak_ddm_counts(j,i);
        peak_signal_watts_LHCP1 = peak_ddm_watts(j,i);
        
        peak_power_counts_RHCP1 = peak_ddm_counts(j+J/2,i);
        peak_signal_watts_RHCP1 = peak_ddm_watts(j+J/2,i);

        if ~isnan(peak_power_counts_LHCP1)

            % derive SNR
            snr_LHCP1 = peak_signal_watts_LHCP1/avg_noise_watts_LHCP1;
            snr_db_LHCP1 = pow2db(snr_LHCP1);

            snr_RHCP1 = peak_signal_watts_RHCP1/avg_noise_watts_RHCP1;
            snr_db_RHCP1 = pow2db(snr_RHCP1);

            % derive instrument gain
            peak_signal_counts_LHCP1 = peak_power_counts_LHCP1-avg_noise_counts_LHCP1;
            inst_gain_LHCP1 = peak_signal_counts_LHCP1/peak_signal_watts_LHCP1;

            peak_signal_counts_RHCP1 = peak_power_counts_RHCP1-avg_noise_counts_RHCP1;
            inst_gain_RHCP1 = peak_signal_counts_RHCP1/peak_signal_watts_RHCP1;

            % save variables
            noise_floor_counts(j,i) = avg_noise_counts_LHCP1;
            noise_floor_counts(j+J/2,i) = avg_noise_counts_RHCP1;

            noise_floor(j,i) = avg_noise_watts_LHCP1;
            noise_floor(j+J/2,i) = avg_noise_watts_RHCP1;

            inst_gain(j,i) = inst_gain_LHCP1;
            inst_gain(j+J/2,i) = inst_gain_RHCP1;

            snr_db(j,i) = snr_db_LHCP1;
            snr_db(j+J/2,i) = snr_db_RHCP1;

        end
    end
end
%}
% save outputs to L1 structure
L1_postCal.raw_counts = ddm_power_counts;
L1_postCal.l1a_power_ddm = power_analog;
L1_postCal.zenith_sig_i2q2 = zenith_i2q2;

%L1_postCal.ddm_noise_floor = noise_floor;
%L1_postCal.ddm_snr = snr_db;

%L1_postCal.inst_gain = inst_gain;
L1_postCal.ddm_ant = ddm_ant;

% Part 3 ends

% Part 4A: SP solver and geometries

% initialise variables
sx_pos_x = zeros(J,I)+invalid;
sx_pos_y = zeros(J,I)+invalid;
sx_pos_z = zeros(J,I)+invalid;

sx_lat = zeros(J,I)+invalid;
sx_lon = zeros(J,I)+invalid;
sx_alt = zeros(J,I)+invalid;

sx_vel_x = zeros(J,I)+invalid;
sx_vel_y = zeros(J,I)+invalid;
sx_vel_z = zeros(J,I)+invalid;

sx_inc_angle = zeros(J,I)+invalid;
sx_d_snell_angle = zeros(J,I)+invalid;
dist_to_coast_km = zeros(J,I)+invalid;
surface_type = zeros(J,I)+invalid;

LOS_flag = zeros(J,I)+invalid;

tx_to_sp_range = zeros(J,I)+invalid;
rx_to_sp_range = zeros(J,I)+invalid;

gps_boresight = zeros(J,I)+invalid;

sx_theta_body = zeros(J,I)+invalid;
sx_az_body = zeros(J,I)+invalid;

sx_theta_enu = zeros(J,I)+invalid;
sx_az_enu = zeros(J,I)+invalid;

gps_tx_power_db_w = zeros(J,I)+invalid;
gps_ant_gain_db_i = zeros(J,I)+invalid;
static_gps_eirp = zeros(J,I)+invalid;

sx_rx_gain_copol = zeros(J,I)+invalid;
sx_rx_gain_xpol = zeros(J,I)+invalid;

for i = 1:I

    % retrieve rx positions, velocities and attitdues
    rx_pos_xyz1 = rx_pos_xyz(i,:);      rx1.rx_pos_xyz = rx_pos_xyz1;
    rx_vel_xyz1 = rx_vel_xyz(i,:);      rx1.rx_vel_xyz = rx_vel_xyz1;
    rx_attitude1 = rx_attitude(i,:);    rx1.rx_attitude = rx_attitude1;    

    % variables are solved only for LHCP channels
    % RHCP channels share the same vales except RX gain solved for each
    % channel
    for j = 1:J/2

        % retrieve tx positions and velocities
        tx_pos_xyz1 = [tx_pos_x(j,i) tx_pos_y(j,i) tx_pos_z(j,i)];
        tx_vel_xyz1 = [tx_vel_x(j,i) tx_vel_y(j,i) tx_vel_z(j,i)];
        tx1.tx_pos_xyz = tx_pos_xyz1;
        tx1.tx_vel_xyz = tx_vel_xyz1;

        %trans_id1 = prn_code(j,i);
        sv_num1 = sv_num(j,i);          tx1.sv_num = sv_num1;

        %ddm_ant1 = ddm_ant(j,i);

        % only process these with valid TX positions
        if ~isnan(tx_pos_x(j,i))

            % Part 4.1: SP solver
            % derive SP positions, angle of incidence and distance
            % to coast
            [sx_pos_xyz1,inc_angle_deg1,d_snell_deg1,dist_to_coast_km1,LOS_flag1] = sp_solver(tx_pos_xyz1,rx_pos_xyz1, ...
                dem,dtu10,landmask_nz);

            LOS_flag(j,i) = LOS_flag1;

            % only process samples with valid sx positions, i.e., LOS = 1
            if LOS_flag1 == 1

                sx_pos_lla1 = ecef2lla(sx_pos_xyz1);            % <lat,lon,alt> of the specular reflection
                surface_type1 = get_surf_type2(sx_pos_xyz1,landmask_nz,water_mask,lcv_mask);    % algorithm version 1.11

                % derive sx velocity
                dt = 1;                                         % time step in second
                tx_pos_xyz_dt = tx_pos_xyz1+dt*tx_vel_xyz1;
                rx_pos_xyz_dt = rx_pos_xyz1+dt*rx_vel_xyz1;
                [sx_pos_xyz_dt,~,~,~,~] = sp_solver(tx_pos_xyz_dt,rx_pos_xyz_dt,dem,dtu10,landmask_nz);

                sx_vel_xyz1 = sx_pos_xyz_dt-sx_pos_xyz1;            

                % save sx values to variables
                sx_pos_x(j,i) = sx_pos_xyz1(1);
                sx_pos_y(j,i) = sx_pos_xyz1(2);
                sx_pos_z(j,i) = sx_pos_xyz1(3);

                sx_lat(j,i) = sx_pos_lla1(1);
                sx_lon(j,i) = sx_pos_lla1(2);
                sx_alt(j,i) = sx_pos_lla1(3);

                sx_vel_x(j,i) = sx_vel_xyz1(1);
                sx_vel_y(j,i) = sx_vel_xyz1(2);
                sx_vel_z(j,i) = sx_vel_xyz1(3);
                surface_type(j,i) = surface_type1;

                sx_inc_angle(j,i) = inc_angle_deg1;
                sx_d_snell_angle(j,i) = d_snell_deg1;
                dist_to_coast_km(j,i) = dist_to_coast_km1;            

                % Part 4.2: SP-related variables - 1
                % this part derives tx/rx gains, ranges and other related
                % variables
                
                % derive SP related geo-parameters, including angles
                % in various frames, ranges and antenna gain/GPS EIRP
                [sx_angle_body1,sx_angle_enu1,sx_angle_ant1,theta_gps1,ranges1,gps_rad1] = spRelated(tx1,rx1, ...
                    sx_pos_xyz1,SV_eirp_LUT);

                % get active antenna gain for LHCP and RHCP channels
                sx_rx_gain_LHCP1 = get_sx_rx_gain(sx_angle_ant1,LHCP_pattern);
                sx_rx_gain_RHCP1 = get_sx_rx_gain(sx_angle_ant1,RHCP_pattern);

                % save to variables
                sx_theta_body(j,i) = sx_angle_body1(1);
                sx_az_body(j,i) = sx_angle_body1(2);

                sx_theta_enu(j,i) = sx_angle_enu1(1);
                sx_az_enu(j,i) = sx_angle_enu1(2);

                gps_boresight(j,i) = theta_gps1;

                tx_to_sp_range(j,i) = ranges1(1);
                rx_to_sp_range(j,i) = ranges1(2);

                gps_tx_power_db_w(j,i) = gps_rad1(1);
                gps_ant_gain_db_i(j,i) = gps_rad1(2);
                static_gps_eirp(j,i) = gps_rad1(3);

                sx_rx_gain_copol(j,i) = sx_rx_gain_LHCP1(1);      % LHCP channel LHCP rx gain
                sx_rx_gain_copol(j+J/2,i) = sx_rx_gain_RHCP1(2);  % RHCP channel RHCP rx gain

                % algorithm version 1.11 changes
                sx_rx_gain_xpol(j,i) = sx_rx_gain_LHCP1(2);       % LHCP channel RHCP rx gain
                sx_rx_gain_xpol(j+J/2,i) = sx_rx_gain_RHCP1(1);   % RHCP channel LHCP rx gain
                                              
            end
        end

    end
end

% expand to RHCP channels
sx_pos_x(J/2+1:J,:) = sx_pos_x(1:J/2,:);
sx_pos_y(J/2+1:J,:) = sx_pos_y(1:J/2,:);
sx_pos_z(J/2+1:J,:) = sx_pos_z(1:J/2,:);

sx_lat(J/2+1:J,:) = sx_lat(1:J/2,:);
sx_lon(J/2+1:J,:) = sx_lon(1:J/2,:);
sx_alt(J/2+1:J,:) = sx_alt(1:J/2,:);

sx_vel_x(J/2+1:J,:) = sx_vel_x(1:J/2,:);
sx_vel_y(J/2+1:J,:) = sx_vel_y(1:J/2,:);
sx_vel_z(J/2+1:J,:) = sx_vel_z(1:J/2,:);

surface_type(J/2+1:J,:) = surface_type(1:J/2,:);
dist_to_coast_km(J/2+1:J,:) = dist_to_coast_km(1:J/2,:);
LOS_flag(J/2+1:J,:) = LOS_flag(1:J/2,:);

rx_to_sp_range(J/2+1:J,:) = rx_to_sp_range(1:J/2,:);
tx_to_sp_range(J/2+1:J,:) = tx_to_sp_range(1:J/2,:);

sx_inc_angle(J/2+1:J,:) = sx_inc_angle(1:J/2,:);
sx_d_snell_angle(J/2+1:J,:) = sx_d_snell_angle(1:J/2,:);

sx_theta_body(J/2+1:J,:) = sx_theta_body(1:J/2,:);
sx_az_body(J/2+1:J,:) = sx_az_body(1:J/2,:);

sx_theta_enu(J/2+1:J,:) = sx_theta_enu(1:J/2,:);
sx_az_enu(J/2+1:J,:) = sx_az_enu(1:J/2,:);

gps_boresight(J/2+1:J,:) = gps_boresight(1:J/2,:);

static_gps_eirp(J/2+1:J,:) = static_gps_eirp(1:J/2,:);

gps_tx_power_db_w(J/2+1:J,:) = gps_tx_power_db_w(1:J/2,:);
gps_ant_gain_db_i(J/2+1:J,:) = gps_ant_gain_db_i(1:J/2,:);

% derive noise floor, SNR and instrument gain
% new change in ver 1.21
snr_LHCP = zeros(J/2,I)+invalid;
snr_RHCP = zeros(J/2,I)+invalid;

snr_flag_LHCP = zeros(J/2,I)+invalid;
snr_flag_RHCP = zeros(J/2,I)+invalid;

inst_gain_LHCP = zeros(J/2,I)+invalid;
inst_gain_RHCP = zeros(J/2,I)+invalid;

noise_floor_LHCP_all = ddm_noise_counts(1:J/2,:);
noise_floor_RHCP_all = ddm_noise_counts(J/2+1:J,:);

peak_delay_row_LHCP = peak_delay_bin(1:J/2,:);
peak_delay_row_RHCP = peak_delay_bin(J/2+1:J,:);

% evaluate new noise floor
dist_to_coast_LHCP = dist_to_coast_km(1:10,:);
ocean_idx = find(dist_to_coast_LHCP<-5);

if ~isempty(ocean_idx)          % with ocean DDMs
    noise_floor_LHCP_ocean = noise_floor_LHCP_all(ocean_idx);
    noise_floor_LHCP = mean(noise_floor_LHCP_ocean(noise_floor_LHCP_ocean>0));

    noise_floor_RHCP_ocean = noise_floor_RHCP_all(ocean_idx);
    noise_floor_RHCP = mean(noise_floor_RHCP_ocean(noise_floor_RHCP_ocean>0));

elseif isempty(ocean_idx)       % without ocean DDMs
    peak_row_idx_LHCP = peak_delay_row_LHCP<30;
    noise_floor_LHCP_land = noise_floor_LHCP_all(peak_row_idx_LHCP);
    noise_floor_LHCP = mean(noise_floor_LHCP_land);

    peak_row_idx_RHCP = peak_delay_row_RHCP<30;
    noise_floor_RHCP_land = noise_floor_RHCP_all(peak_row_idx_RHCP);
    noise_floor_RHCP = mean(noise_floor_RHCP_land);

end

for j = 1:J/2
    for i = 1:I

        peak_counts_LHCP1 = peak_ddm_counts(j,i);
        peak_counts_RHCP1 = peak_ddm_counts(j+J/2,i);

        peak_signal_watts_LHCP1 = peak_ddm_watts(j,i);
        peak_signal_watts_RHCP1 = peak_ddm_watts(j+J/2,i);

        if ~isnan(peak_counts_LHCP1)
            peak_signal_LHCP1 = peak_counts_LHCP1-noise_floor_LHCP;

            if peak_signal_LHCP1>0           
                snr_LHCP1 = peak_signal_LHCP1/noise_floor_LHCP;               
                snr_flag_LHCP1 = 0;
            
            elseif peak_signal_LHCP1<=0
                snr_LHCP1 = nan;
                snr_flag_LHCP1 = 1;
                
            end
            
            snr_LHCP(j,i) = pow2db(snr_LHCP1);
            snr_flag_LHCP(j,i) = snr_flag_LHCP1;

        end

        if ~isnan(peak_counts_RHCP1)
            peak_signal_RHCP1 = peak_counts_RHCP1-noise_floor_RHCP;
            
            if peak_signal_RHCP1>0
                snr_RHCP1 = peak_signal_RHCP1/noise_floor_RHCP; 
                snr_flag_RHCP1 = 0;
          
            elseif peak_signal_RHCP1<=0
                snr_RHCP1 = nan;
                snr_flag_RHCP1 = 1;
            
            end

            snr_RHCP(j,i) = pow2db(snr_RHCP1);
            snr_flag_RHCP(j,i) = snr_flag_RHCP1;

        end

        inst_gain_LHCP1 = peak_counts_LHCP1/peak_signal_watts_LHCP1;
        inst_gain_RHCP1 = peak_counts_RHCP1/peak_signal_watts_RHCP1;

        inst_gain_LHCP(j,i) = inst_gain_LHCP1;
        inst_gain_RHCP(j,i) = inst_gain_RHCP1;

    end
end

ddm_snr = [snr_LHCP;snr_RHCP];
ddm_noise_floor = [repmat(noise_floor_LHCP,[10,I]);repmat(noise_floor_RHCP,[10,I])];
ddm_snr_flag = [snr_flag_LHCP;snr_flag_RHCP];

inst_gain = [inst_gain_LHCP;inst_gain_RHCP];

% save variables
L1_postCal.sp_pos_x = sx_pos_x;
L1_postCal.sp_pos_y = sx_pos_y;
L1_postCal.sp_pos_z = sx_pos_z;

L1_postCal.sp_lat = sx_lat;
L1_postCal.sp_lon = sx_lon;
L1_postCal.sp_alt = sx_alt;

L1_postCal.sp_vel_x = sx_vel_x;
L1_postCal.sp_vel_y = sx_vel_y;
L1_postCal.sp_vel_z = sx_vel_z;

L1_postCal.sp_surface_type = surface_type;
L1_postCal.sp_dist_to_coast_km = dist_to_coast_km;
L1_postCal.LOS_flag = LOS_flag;

L1_postCal.rx_to_sp_range = rx_to_sp_range;
L1_postCal.tx_to_sp_range = tx_to_sp_range;

L1_postCal.sp_inc_angle = sx_inc_angle;
L1_postCal.sp_d_snell_angle = sx_d_snell_angle;

L1_postCal.sp_theta_body = sx_theta_body;
L1_postCal.sp_az_body = sx_az_body;
L1_postCal.sp_theta_enu = sx_theta_enu;
L1_postCal.sp_az_enu = sx_az_enu;

L1_postCal.sp_rx_gain_copol = sx_rx_gain_copol;
L1_postCal.sp_rx_gain_xpol = sx_rx_gain_xpol;       % algorithm version 1.11

L1_postCal.gps_off_boresight_angle_deg = gps_boresight;

L1_postCal.static_gps_eirp = static_gps_eirp;
L1_postCal.gps_tx_power_db_w = gps_tx_power_db_w;
L1_postCal.gps_ant_gain_db_i = gps_ant_gain_db_i;

L1_postCal.ddm_noise_floor = ddm_noise_floor;
L1_postCal.ddm_snr = ddm_snr;
L1_postCal.ddm_snr_flag = ddm_snr_flag;

L1_postCal.inst_gain = inst_gain;

% Part 4B: BRCS/NBRCS, reflectivity, coherent status and fresnel zone
%clc

% initialise variables
brcs_ddm_peak_bin_delay_row = zeros(J,I)+invalid;
brcs_ddm_peak_bin_dopp_col = zeros(J,I)+invalid;

brcs_ddm_sp_bin_delay_row = zeros(J,I)+invalid;
brcs_ddm_sp_bin_dopp_col = zeros(J,I)+invalid;

sp_delay_error = zeros(J,I)+invalid;
sp_dopp_error = zeros(J,I)+invalid;

confidence_flag = zeros(J,I)+invalid;

zenith_code_phase = zeros(J,I)+invalid;

brcs = zeros(5,40,J,I)+invalid;
A_eff = zeros(5,40,J,I)+invalid;
A_eff_all = zeros(9,79,J,I)+invalid;            % debug only

norm_refl_waveform = zeros(1,40,J,I)+invalid;

nbrcs_scatter_area_v1 = zeros(J,I)+invalid;
ddm_nbrcs_v1 = zeros(J,I)+invalid;

nbrcs_scatter_area_v2 = zeros(J,I)+invalid;
ddm_nbrcs_v2 = zeros(J,I)+invalid;

%les_scatter_area = zeros(J,I)+invalid;
%ddm_les = zeros(J,I)+invalid;

%tes_scatter_area = zeros(J,I)+invalid;
%ddm_tes = zeros(J,I)+invalid;

surface_reflectivity = zeros(5,40,J,I)+invalid;
surface_reflectivity_peak = zeros(J,I)+invalid;

fresnel_coeff = zeros(J,I)+invalid;
fresnel_minor = zeros(J,I)+invalid;
fresnel_major = zeros(J,I)+invalid;
fresnel_orientation = zeros(J,I)+invalid;

coherency_ratio = zeros(J,I)+invalid;
coherency_state = zeros(J,I)+invalid;

% derive amb-function (chi2) to be used in computing A_eff
chi2 = get_chi2(40,5);

% derive floating SP bin location and effective scattering area A_eff
for i = 1:I

    % retrieve rx positions and velocities
    rx_pos_xyz1 = rx_pos_xyz(i,:);          rx1.rx_pos_xyz = rx_pos_xyz1;
    rx_vel_xyz1 = rx_vel_xyz(i,:);          rx1.rx_vel_xyz = rx_vel_xyz1;
    rx_clk_drift1 = rx_clk_drift_mps(i,:);  rx1.rx_clk_drift = rx_clk_drift1;

    for j = 1:J/2

        % retrieve tx positions and velocities
        tx_pos_xyz1 = [tx_pos_x(j,i) tx_pos_y(j,i) tx_pos_z(j,i)];
        tx_vel_xyz1 = [tx_vel_x(j,i) tx_vel_y(j,i) tx_vel_z(j,i)];
        tx1.tx_pos_xyz = tx_pos_xyz1;
        tx1.tx_vel_xyz = tx_vel_xyz1;

        % retrieve sx-related parameters
        sx_pos_xyz1 = [sx_pos_x(j,i) sx_pos_y(j,i) sx_pos_z(j,i)];
        sx_pos_lla1 = ecef2lla(sx_pos_xyz1);

        %sx_inc_angle1 = sx_inc_angle(j,i);
        sx_d_snell_deg1 = sx_d_snell_angle(j,i);
        dist_to_coast1 = dist_to_coast_km(j,i);
        
        sx1.sx_pos_xyz = sx_pos_xyz1;
        sx1.sx_d_snell = sx_d_snell_deg1;
        sx1.dist_to_coast = dist_to_coast1;

        % retrieve ddm-related variables
        raw_counts1 = raw_counts(:,:,j,i);          ddm1.raw_counts = raw_counts1;        
        add_range_to_sp1 = add_range_to_sp(j,i);    ddm1.add_range_to_sp = add_range_to_sp1;
        snr_db1 = snr_db(j,i);                      ddm1.snr_db = snr_db1;

        delay_center_chips1 = delay_center_chips(j,i);
        ddm1.delay_center_chips = delay_center_chips1;

        doppler_center_hz1 = doppler_center_hz(j,i);
        ddm1.doppler_center_hz = doppler_center_hz1;

        T_coh1 = coherent_duration(i);
        ddm1.T_coh = T_coh1;

        ddm1.delay_resolution = 0.25;   ddm1.num_delay_bins = 40;   ddm1.delay_center_bin = 20;
        ddm1.doppler_resolution = 500;  ddm1.num_doppler_bins = 5;  ddm1.doppler_center_bin = 2;

        if ~isnan(sx_pos_x(j,i)) && (sum(raw_counts1,'all')~=0)
            
            % Part 4.3: SP-related variables - 2
            % this part derives confidence and floating bin locations of SP
            raw_counts_max = max(raw_counts1,[],'all');
            [peak_doppler_bin1,peak_delay_bin1] = find(raw_counts1 == raw_counts_max,1);

            [specular_bin1,zenith_code_phase1,confidence_flag1] = get_specular_bin(tx1,rx1,sx1,ddm1);

            sx1.sx_delay_bin = specular_bin1(1);
            sx1.sx_doppler_bin = specular_bin1(2);

            % Part 4.4a: Effective scattering area
            L = 18030; grid_res = 30;       % L may need to be updated in the future
            local_dem1 = get_local_dem(sx_pos_lla1,L,grid_res,dem,dtu10,dist_to_coast1);

            [A_eff1,A_eff_all1] = get_ddm_Aeff(tx1,rx1,sx1,local_dem1,phy_ele_size,chi2);

            % save to variables
            brcs_ddm_peak_bin_delay_row(j,i) = peak_delay_bin1-1;   % minus 1 for 0-based indces
            brcs_ddm_peak_bin_dopp_col(j,i) = peak_doppler_bin1-1;

            brcs_ddm_sp_bin_delay_row(j,i) = specular_bin1(1);
            brcs_ddm_sp_bin_dopp_col(j,i) = specular_bin1(2);
            sp_delay_error(j,i) = specular_bin1(3);
            sp_dopp_error(j,i) = specular_bin1(4);
            
            zenith_code_phase(j,i) = zenith_code_phase1;

            confidence_flag(j,i) = confidence_flag1;

            A_eff(:,:,j,i) = A_eff1;
            A_eff_all(:,:,j,i) = A_eff_all1;

        end

    end
end

% extend to RHCP channels
brcs_ddm_peak_bin_delay_row(J/2+1:J,:) = brcs_ddm_peak_bin_delay_row(1:J/2,:);
brcs_ddm_peak_bin_dopp_col(J/2+1:J,:) = brcs_ddm_peak_bin_dopp_col(1:J/2,:);

brcs_ddm_sp_bin_delay_row(J/2+1:J,:) = brcs_ddm_sp_bin_delay_row(1:J/2,:);
brcs_ddm_sp_bin_dopp_col(J/2+1:J,:) = brcs_ddm_sp_bin_dopp_col(1:J/2,:);
sp_delay_error(J/2+1:J,:) = sp_delay_error(1:J/2,:);
sp_dopp_error(J/2+1:J,:) = sp_dopp_error(1:J/2,:);
            
zenith_code_phase(J/2+1:J,:) = zenith_code_phase(1:J/2,:);

confidence_flag(J/2+1:J,:) = confidence_flag(1:J/2,:);

A_eff(:,:,J/2+1:J,:) = A_eff(:,:,1:J/2,:);
A_eff_all(:,:,J/2+1:J,:) = A_eff_all(:,:,1:J/2,:);

% save variables
L1_postCal.brcs_ddm_peak_bin_delay_row = brcs_ddm_peak_bin_delay_row;
L1_postCal.brcs_ddm_peak_bin_dopp_col = brcs_ddm_peak_bin_dopp_col;

L1_postCal.brcs_ddm_sp_bin_delay_row = brcs_ddm_sp_bin_delay_row;
L1_postCal.brcs_ddm_sp_bin_dopp_col = brcs_ddm_sp_bin_dopp_col;

L1_postCal.sp_delay_error = sp_delay_error;
L1_postCal.sp_dopp_error = sp_dopp_error;  
L1_postCal.sp_ngrx_delay_correction = sp_delay_error;
L1_postCal.sp_ngrx_dopp_correction = sp_dopp_error;

L1_postCal.zenith_code_phase = zenith_code_phase;

L1_postCal.confidence_flag = confidence_flag;

L1_postCal.eff_scatter = A_eff;
L1_postCal.A_eff_all = A_eff_all;

% derive brcs, nbrcs, and other parameters
for i = 1:I
    for j = 1:J

        % variables for deriving BRCS and reflectivity
        tx_pos_xyz1 = [tx_pos_x(j,i) tx_pos_y(j,i) tx_pos_z(j,i)];
        rx_pos_xyz1 = rx_pos_xyz(i,:);
        sx_pos_xyz1 = [sx_pos_x(j,i) sx_pos_y(j,i) sx_pos_z(j,i)];

        inc_angle1 = sx_inc_angle(j,i);
        dist_to_coast1 = dist_to_coast_km(j,i);

        eirp_watt1 = static_gps_eirp(j,i);
        rx_gain_db_i1 = sx_rx_gain_copol(j,i);
        TSx1 = tx_to_sp_range(j,i);
        RSx1 = rx_to_sp_range(j,i);
        
        ddm_ant1 = ddm_ant(j,i);
        
        % retrieve ddm-related variables
        raw_counts1 = ddm_power_counts(:,:,j,i); 
        power_analog1 = power_analog(:,:,j,i);
        snr_db1 = snr_db(j,i);       
        
        if ~isnan(ddm_ant1) && ~isnan(sx_pos_x(j,i)) && (sum(raw_counts1,'all')~=0)

            % compensate cable loss
            if ddm_ant1 == 2
                cable_loss_db = 0.6600;         % LHCP cable loss

            elseif ddm_ant1 == 3
                cable_loss_db = 0.5840;         % RHCP cable loss

            end

            cable_loss = db2pow(cable_loss_db);
            power_analog_cable_loss1 = power_analog1*cable_loss;
            
            % Part 4.4b: brcs, nbrcs
            brcs1 = ddm_brcs(power_analog_cable_loss1,eirp_watt1,rx_gain_db_i1,TSx1,RSx1);

            A_eff1 = A_eff(:,:,j,i);
            sx_bin1(1) = brcs_ddm_sp_bin_delay_row(j,i);
            sx_bin1(2) = brcs_ddm_sp_bin_dopp_col(j,i);
            
            % the below computes two versions of NBRCS
            % version 1: smaller area, version 2: larger area
            [nbrcs_v1_1,nbrcs_scatter_v1_1] = get_ddm_nbrcs2(brcs1,A_eff1,sx_bin1,1);
            [nbrcs_v2_1,nbrcs_scatter_v2_1] = get_ddm_nbrcs2(brcs1,A_eff1,sx_bin1,2);

            % Part 4.5: reflectivity and peak reflectivity
            [refl1,refl_peak1] = ddm_refl(power_analog_cable_loss1,eirp_watt1,rx_gain_db_i1,TSx1,RSx1);
           
            % Part 4.6: Fresnel coefficient and dimensions
            [fresnel_coeff1,fresnel_axis1,fresnel_orientation1] = get_fresnel(tx_pos_xyz1, ...
                rx_pos_xyz1,sx_pos_xyz1,dist_to_coast1,inc_angle1,ddm_ant1);

            % Part 4.7: coherent status
            [CR1,CS1] = coh_det(raw_counts1,snr_db1);

            % normalised reflected waveform
            refl_waveform1 = sum(refl1,1);
            norm_refl_waveform1 = refl_waveform1/refl_peak1;
            
            % save to variables
            brcs(:,:,j,i) = brcs1;

            nbrcs_scatter_area_v1(j,i) = nbrcs_scatter_v1_1;
            ddm_nbrcs_v1(j,i) = nbrcs_v1_1;

            nbrcs_scatter_area_v2(j,i) = nbrcs_scatter_v2_1;
            ddm_nbrcs_v2(j,i) = nbrcs_v2_1;
            
            %nbrcs_scatter_area(j,i) = nbrcs1.nbrcs_scatter;
            %ddm_nbrcs(j,i) = nbrcs1.nbrcs_value;

            %les_scatter_area(j,i) = LES1.LES_scatter;
            %ddm_les(j,i) = LES1.LES_slope;

            %tes_scatter_area(j,i) = TES1.TES_scatter;
            %ddm_tes(j,i) = TES1.TES_slope;

            surface_reflectivity(:,:,j,i) = refl1;
            surface_reflectivity_peak(j,i) = refl_peak1;

            fresnel_coeff(j,i) = fresnel_coeff1;
            fresnel_major(j,i) = fresnel_axis1(1);
            fresnel_minor(j,i) = fresnel_axis1(2);
            fresnel_orientation(j,i) = fresnel_orientation1;

            coherency_ratio(j,i) = CR1;
            coherency_state(j,i) = CS1;

            norm_refl_waveform(:,:,j,i) = norm_refl_waveform1;
            
        end
        
    end
end

L1_postCal.brcs = brcs;

L1_postCal.nbrcs_scatter_area_v1 = nbrcs_scatter_area_v1;
L1_postCal.ddm_nbrcs_v1 = ddm_nbrcs_v1;

L1_postCal.nbrcs_scatter_area_v2 = nbrcs_scatter_area_v2;
L1_postCal.ddm_nbrcs_v2 = ddm_nbrcs_v2;

%L1_postCal.les_scatter_area = les_scatter_area;
%L1_postCal.ddm_les = ddm_les;

%L1_postCal.tes_scatter_area = tes_scatter_area;
%L1_postCal.ddm_tes = ddm_tes;

L1_postCal.surface_reflectivity = surface_reflectivity;
L1_postCal.surface_reflectivity_peak = surface_reflectivity_peak;

L1_postCal.fresnel_coeff = fresnel_coeff;
L1_postCal.fresnel_major = fresnel_major;
L1_postCal.fresnel_minor = fresnel_minor;
L1_postCal.fresnel_orientation = fresnel_orientation;

L1_postCal.coherency_ratio = coherency_ratio;
L1_postCal.coherency_state = coherency_state;

L1_postCal.norm_refl_waveform = norm_refl_waveform;

% Cross Pol
nbrcs_cross_pol_v1 = zeros(J,I)+invalid;
nbrcs_cross_pol_v2 = zeros(J,I)+invalid;

for i = 1:I
    for j = 1:J/2

        nbrcs_LHCP_v1 = ddm_nbrcs_v1(j,i);
        nbrcs_RHCP_v1 = ddm_nbrcs_v1(j+J/2,i);

        nbrcs_LHCP_v2 = ddm_nbrcs_v2(j,i);
        nbrcs_RHCP_v2 = ddm_nbrcs_v2(j+J/2,i);

        CP1 = nbrcs_LHCP_v1/nbrcs_RHCP_v1;
        CP_db1 = pow2db(CP1);

        CP2 = nbrcs_LHCP_v2/nbrcs_RHCP_v2;
        CP_db2 = pow2db(CP2);

        nbrcs_cross_pol_v1(j,i) = CP_db1;
        nbrcs_cross_pol_v2(j,i) = CP_db2;

    end
end

nbrcs_cross_pol_v1(11:20,:) = -1*nbrcs_cross_pol_v1(1:10,:);
nbrcs_cross_pol_v2(11:20,:) = -1*nbrcs_cross_pol_v2(1:10,:);

L1_postCal.nbrcs_cross_pol_v1 = nbrcs_cross_pol_v1;
L1_postCal.nbrcs_cross_pol_v2 = nbrcs_cross_pol_v2;
L1_postCal.lna_noise_figure = zeros(J,I)+3;         % LNA noise figure is 3 dB according to the specification

% Quality Flags
quality_flags1 = zeros(J,I)+invalid;

for i = 1:I
    for j = 1:J

        quality_flag1_1 = zeros(1,24);

        % flag 2, 3 and 23
        rx_roll1 = rx_roll(i);
        rx_pitch1 = rx_pitch(i);
        rx_yaw1 = rx_yaw(i);
    
        if (rx_roll1 >= 30) || (rx_pitch1 >= 10) || (rx_yaw1 >= 5)
            quality_flag1_1(3) = 1;
        else
            quality_flag1_1(2) = 1;
        end

        if rx_roll1 > 1
            quality_flag1_1(23) = 1;
        end

        % flag 4
        quality_flag1_1(4) = 0;

        % flag 5 and 6
        trans_id1 = transmitter_id(j,i);
        if trans_id1 == 0
            quality_flag1_1(5) = 1;
        end

        if trans_id1 == 28
            quality_flag1_1(6) = 1;
        end

        % flag 7 and 10
        snr_db1 = snr_db(j,i);

        if i > 1
            snr_db2 = snr_db(j,i-1);
            diff1 = (db2pow(snr_db1)-db2pow(snr_db2))/db2pow(snr_db1);
            diff2 = snr_db1-snr_db2;

            if abs(diff1) > 0.1
                quality_flag1_1(7) = 1;
            end

            if abs(diff2) > 0.24
                quality_flag1_1(10) = 1;
            end
        end

        % flag 8 and 9
        dist_to_coast1 = dist_to_coast_km(j,i);

        if dist_to_coast1 > 0
            quality_flag1_1(8) = 1;
        end

        if dist_to_coast1 > -25
            quality_flag1_1(9) = 1;
        end

        % flag 11
        ant_temp1 = ant_temp_nadir(i);
        if i > 1
            ant_temp2 = ant_temp_nadir(i-1);
            rate = (ant_temp2-ant_temp1)*60;

            if rate > 1
                quality_flag1_1(11) = 1;
            end

        end

        % flag 12
        zenith_code_phase1 = zenith_code_phase(j,i);
        signal_code_phase1 = delay_correction(meter2chips(add_range_to_sp(j,i)),1023);
        diff1 = zenith_code_phase1-signal_code_phase1;
        if diff1 >= 10
            quality_flag1_1(12) = 1;
        end

        % flag 15 and 16
        sp_delay_row = brcs_ddm_sp_bin_delay_row(j,i);
        sp_dopp_col = brcs_ddm_sp_bin_dopp_col(j,i);

        if (sp_delay_row<15)||(sp_delay_row>35)
            quality_flag1_1(15) = 1;
        end

        if (sp_dopp_col<2)||(sp_dopp_col>4)
            quality_flag1_1(16) = 1;
        end

        % flag 17
        if (floor(sp_delay_row) < 38) && (floor(sp_delay_row) > 0) && ...
                (floor(sp_dopp_col) < 5) && (floor(sp_dopp_col) > 1)
            brcs_ddma = brcs(floor(sp_dopp_col)-1:floor(sp_dopp_col)+1, ...
                            floor(sp_delay_row):floor(sp_dopp_col)+3);
            det = find(brcs_ddma<0, 1);
            if ~isempty(det)
                quality_flag1_1(17) = 1;
            end
        end

        % flag 18
        tx_pos_x1 = tx_pos_x(j,i);
        prn_code1 = prn_code(j,i);
        if (tx_pos_x1 == 0) && (~isnan(prn_code1))
            quality_flag1_1(18) = 1;
        end

        % flag 19
        sx_pos_x1 = sx_pos_x(j,i);
        if (sx_pos_x1 == invalid) && (~isnan(prn_code1))
            quality_flag1_1(19) = 1;
        end

        % flag 20
        rx_gain1 = sx_rx_gain_copol(j,i);
        if (rx_gain1 == invalid) && (~isnan(prn_code1))
            quality_flag1_1(20) = 1;
        end

        quality_flag1_1(21) = 1;

        % flag 22
        rx_alt = rx_pos_lla(i,3);
        if rx_alt > 15000
            quality_flag1_1(22) = 1;
        end

        % flag 24
        prn1 = prn_code(j,i);
        if prn1 == 28
            quality_flag1_1(24) = 1;
        end

        % flag 1
        if (quality_flag1_1(3) == 1 || ...
            quality_flag1_1(4) == 1 || ...
            quality_flag1_1(5) == 1 || ...
            quality_flag1_1(6) == 1 || ...
            quality_flag1_1(7) == 1 || ...
            quality_flag1_1(10) == 1 || ...
            quality_flag1_1(11) == 1 || ...
            quality_flag1_1(12) == 1 || ...
            quality_flag1_1(13) == 1 || ...
            quality_flag1_1(14) == 1 || ...
            quality_flag1_1(15) == 1 || ...
            quality_flag1_1(16) == 1 || ...
            quality_flag1_1(17) == 1 || ...
            quality_flag1_1(18) == 1 || ...
            quality_flag1_1(20) == 1 || ...
            quality_flag1_1(22) == 1 || ...
            quality_flag1_1(23) == 1 || ...
            quality_flag1_1(24) == 1)
        
            quality_flag1_1(1) = 1;

        end

        quality_flags1(j,i) = get_quality_flag(quality_flag1_1);
        
    end
end

L1_postCal.quality_flags1 = quality_flags1;
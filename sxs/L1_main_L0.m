% This script operates the L1b calibration for Rongowai
% Version: Matlab - 1.0
% inital version verified by CYGNSS L1
% Version: Matlab - 1.1
% 1) update according to the Rongowai bench test data
% 2) convert all external source files date files
% 3) the value for all invalid fields is defined as -9999
% 4) L1a changes to using power in watts to avoid error when converting a
% negative power watts to dB

%% Prelaunch 1: Load L0 data
clear
clc

% load L0 netCDF
path1 = '../dat/raw/';
L0_filename = [path1 'Flight_4_2022-10-10.nc'];

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
ddm_number = double(ncread(L0_filename,'/science/ddm/ddm_number'));

% noise bin average
noise_std_dev_rf1 = double(ncread(L0_filename,'/science/ddm/RF1_zenith_RHCP_std_dev'));
noise_std_dev_rf2 = double(ncread(L0_filename,'/science/ddm/RF2_nadir_LHCP_std_dev')); 
noise_std_dev_rf3 = double(ncread(L0_filename,'/science/ddm/RF3_nadir_RHCP_std_dev'));

delay_bin_res = double(ncread(L0_filename,'/science/ddm/delay_bin_res_narrow'));        % delay bin resolution
doppler_bin_res = double(ncread(L0_filename,'/science/ddm/doppler_bin_res_narrow'));    % doppler bin resolution

delay_center_bin = double(ncread(L0_filename,'/science/ddm/ddm_center_delay_bin'));     % ddm center delay bin 
doppler_center_bin = double(ncread(L0_filename,'/science/ddm/ddm_center_doppler_bin')); % ddm center doppler bin

% absolute ddm center delay and doppler
delay_center_chips = double(ncread(L0_filename,'/science/ddm/center_delay_bin_code_phase'));    
doppler_center_hz = double(ncread(L0_filename,'/science/ddm/center_doppler_bin_frequency'));

% number of doppler and delay bins
num_delay_bins = double(ncread(L0_filename,'/science/ddm/num_delay_bins')); 
num_doppler_bins = double(ncread(L0_filename,'/science/ddm/num_doppler_bins'));

% coherent duration and noncoherent integration
coherent_duration = double(ncread(L0_filename,'/science/ddm/L1_E1_coherent_duration'));
non_coherent_integrations = double(ncread(L0_filename,'/science/ddm/L1_E1_non_coherent_integrations'));

% NGRx estimate additional delay path
add_range_to_sp_pvt = double(ncread(L0_filename,'/science/ddm/additional_range_to_SP'));

% antenna temperatures and engineering timestamp
eng_timestamp = double(ncread(L0_filename,'/eng/packet_creation_time'));
zenith_ant_temp_eng = double(ncread(L0_filename,'/eng/zenith_ant_temp'));
nadir_ant_temp_eng = double(ncread(L0_filename,'/eng/nadir_ant_temp'));

% Load L0 data ends

%% Prelaunch 1.5: Filter valid timestampes
clc

c = 299792458;                              % speed of light

% rx-related variables
index1 = ~isnan(pvt_gps_week);

pvt_gps_week = pvt_gps_week(index1);        pvt_gps_sec = pvt_gps_sec(index1);

rx_pos_x_pvt = rx_pos_x_pvt(index1);        rx_pos_y_pvt = rx_pos_y_pvt(index1);
rx_pos_z_pvt = rx_pos_z_pvt(index1);

rx_vel_x_pvt = rx_vel_x_pvt(index1);        rx_vel_y_pvt = rx_vel_y_pvt(index1);
rx_vel_z_pvt = rx_vel_z_pvt(index1);
rx_vel_xyz_pvt = [rx_vel_x_pvt rx_vel_y_pvt rx_vel_z_pvt];

rx_roll_pvt = rx_roll_pvt(index1);          rx_pitch_pvt = rx_pitch_pvt(index1);
rx_yaw_pvt = rx_yaw_pvt(index1);

rx_clk_bias_m_pvt = rx_clk_bias_m_pvt(index1);
rx_clk_drift_mps_pvt = rx_clk_drift_mps_pvt(index1);

% ddm-related variables
index2 = ~isnan(transmitter_id(1,:));

transmitter_id = transmitter_id(:,index2);

first_scale_factor = first_scale_factor(:,index2);
raw_counts = raw_counts(:,:,:,index2); 
zenith_i2q2 = zenith_i2q2(:,index2);

rf_source = rf_source(:,index2);
ddm_number = ddm_number(:,index2);

noise_std_dev_rf1 = noise_std_dev_rf1(index2);
noise_std_dev_rf2 = noise_std_dev_rf2(index2);
noise_std_dev_rf3 = noise_std_dev_rf3(index2);

delay_bin_res = delay_bin_res(index2);      doppler_bin_res = doppler_bin_res(index2);

delay_center_bin = delay_center_bin(index2);
doppler_center_bin = doppler_center_bin(index2);

% absolute ddm center delay and doppler
delay_center_chips = delay_center_chips(:,index2)/c;                    % correction to the chips  
doppler_center_hz = doppler_center_hz(:,index2);

% number of doppler and delay bins
num_delay_bins = num_delay_bins(index2);   num_doppler_bins = num_doppler_bins(index2);

% coherent duration and noncoherent integration
coherent_duration = coherent_duration(index2)/1000;                     % convert to seconds
non_coherent_integrations = non_coherent_integrations(index2)/1000;

% NGRx estimate additional delay path
add_range_to_sp_pvt = add_range_to_sp_pvt(:,index2);

% temperatures from engineering data
index3 = ~isnan(eng_timestamp);

eng_timestamp = eng_timestamp(index3);
nadir_ant_temp_eng = nadir_ant_temp_eng(index3);
zenith_ant_temp_eng = zenith_ant_temp_eng(index3);

%% Prelaunch 2 - define external data paths and filenames
clc

% define IGS orbits filename (*.sp3)
% note this path is defined in C++ format
gps_orbit_filename = '..//dat//orbits//igr22310.sp3';

% load SRTM_30 DEM
dem_path = '../dat/dem/';
dem_file = 'nzsrtm_30_v1.tif';

dem = get_dem([dem_path dem_file]);

% load DTU10 model
dtu_path = '../dat/dtu/';
dtu_filename = 'dtu10_v1.dat';

dtu10 = get_dtu10([dtu_path dtu_filename]);

% load ocean/land (distance to coast) mask
landmask_path = '../dat/cst/';
landmask_filename = 'dist_to_coast_nz_v1.dat';

landmask_nz = get_dist_to_coast_mask([landmask_path landmask_filename]);

% process landcover mask
lcv_path = '../dat/lcv/';
lcv_filename = 'lcv.png';

lcv_mask = imread([lcv_path lcv_filename]);

% process inland water mask
pek_path = '../dat/pek/';

[pek_160E_40S,ref_160E_40S] = readgeoraster([pek_path 'occurrence_160E_40S.tif']);
[pek_170E_30S,ref_170E_30S] = readgeoraster([pek_path 'occurrence_170E_30S.tif']);
[pek_170E_40S,ref_170E_40S] = readgeoraster([pek_path 'occurrence_170E_40S.tif']);

water_mask_160E_40S.pek = pek_160E_40S;
water_mask_160E_40S.ref = ref_160E_40S;

water_mask_170E_30S.pek = pek_170E_30S;
water_mask_170E_30S.ref = ref_170E_30S;

water_mask_170E_40S.pek = pek_170E_40S;
water_mask_170E_40S.ref = ref_170E_40S;

water_mask.water_mask_160E_40S = water_mask_160E_40S;
water_mask.water_mask_170E_30S = water_mask_170E_30S;
water_mask.water_mask_170E_40S = water_mask_170E_40S;

% load PRN-SV and SV-EIRP(static) LUT
gps_path = '../dat/gps/';
SV_PRN_filename = 'PRN_SV_LUT_v1.dat';
SV_eirp_filename = 'GPS_SV_EIRP_Params_v7.dat';

SV_PRN_LUT = readmatrix([gps_path SV_PRN_filename]);
SV_PRN_LUT = SV_PRN_LUT(:,1:2);

SV_eirp_LUT = readmatrix([gps_path SV_eirp_filename]);

% load and process nadir NGRx-GNSS antenna patterns
rng_path = '../dat/rng/';

LHCP_L_filename = 'GNSS_LHCP_L_gain_db_i_v1.dat';
LHCP_R_filename = 'GNSS_LHCP_R_gain_db_i_v1.dat';

RHCP_L_filename = 'GNSS_RHCP_L_gain_db_i_v1.dat';
RHCP_R_filename = 'GNSS_RHCP_R_gain_db_i_v1.dat';

[~,~,LHCP_L_gain_db_i] = get_ant_pattern([rng_path,LHCP_L_filename]);
[~,~,LHCP_R_gain_db_i] = get_ant_pattern([rng_path,LHCP_R_filename]);
LHCP_pattern.LHCP = LHCP_L_gain_db_i;
LHCP_pattern.RHCP = LHCP_R_gain_db_i;

[~,~,RHCP_L_gain_db_i] = get_ant_pattern([rng_path,RHCP_L_filename]);
[~,~,RHCP_R_gain_db_i] = get_ant_pattern([rng_path,RHCP_R_filename]);
RHCP_pattern.LHCP = RHCP_L_gain_db_i;
RHCP_pattern.RHCP = RHCP_R_gain_db_i;

%% Part 1: General processing
% This part derives global constants, timestamps, and all the other
% parameters at ddm timestamps
clc

invalid = -99999999;                                % defines the value to be used for invalid fields

I = length(pvt_gps_sec);                            % total length of samples
J = 20;                                             % maximal NGRx capacity

non_coherent_integrations(isnan(non_coherent_integrations)) = invalid;

% initialise output data array for Part 1 processing
pvt_utc = zeros(I,1)+invalid;       ddm_utc = zeros(I,1)+invalid;
gps_week = zeros(I,1)+invalid;      gps_tow = zeros(I,1)+invalid;
ddm_pvt_bias = zeros(I,1)+invalid;
add_range_to_sp = zeros(J,I)+invalid;
status_flags_one_hz = zeros(I,1)+invalid;

% initialise a structure to save L1 results
L1_postCal = struct;

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
    [year,month,day] = ymd(D_ddm1);
    [hour,mins,secs] = hms(D_ddm1);
    [ddm_gps_week1,ddm_gps_sec1] = utc2gpstime(year,month,day,hour,mins,secs);

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
rx_clk = [rx_clk_bias_m rx_clk_drift_mps];

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

L1_postCal.aircraft_reg = 'ZK-NFA';             % default value
L1_postCal.ddm_source = 2;                      % 1 = GPS signal simulator, 2 = aircraft
L1_postCal.ddm_time_type_selector = 1;          % 1 = middle of DDM sampling period
L1_postCal.delay_resolution = 0.25;             % unit in chips
L1_postCal.dopp_resolution = 500;               % unit in Hz
L1_postCal.dem_source = 'SRTM30';

% write algorithm and LUT versions
L1_postCal.l1_algorithm_version = '1';
L1_postCal.l1_data_version = '1';
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
L1_postCal.pvt_timestamp_utc = pvt_utc;

L1_postCal.ddm_timestamp_gps_week = gps_week;
L1_postCal.ddm_timestamp_gps_sec = gps_tow;
L1_postCal.ddm_timestamp_utc = ddm_utc;

L1_postCal.ddm_pvt_bias = ddm_pvt_bias;

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

%% Part 2: Derive TX related variables
% This part derives TX positions and velocities, maps between PRN and SVN,
% and gets track ID
clc

trans_id_unique = unique(transmitter_id);
trans_id_unique = trans_id_unique(trans_id_unique>0);

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

    for j = 1:J
    
        transmitter_id1 = transmitter_id(j,i);

        % PRN and SVN
        prn1 = transmitter_id1;
        sv_num1 = SV_PRN_LUT(SV_PRN_LUT(:,1)==prn1,2);

        if (~isempty(sv_num1)) && (prn1 ~= 0)
            [tx_pos_xyz1,tx_vel_xyz1,tx_clk_bias1,~] = gps_posvel(prn1,ddm_gps_timestamp1, ...
                gps_orbit_filename);

            tx_pos_x(j,i) = tx_pos_xyz1(1); tx_vel_x(j,i) = tx_vel_xyz1(1);
            tx_pos_y(j,i) = tx_pos_xyz1(2); tx_vel_y(j,i) = tx_vel_xyz1(2);
            tx_pos_z(j,i) = tx_pos_xyz1(3); tx_vel_z(j,i) = tx_vel_xyz1(3);

            tx_clk_bias(j,i) = tx_clk_bias1;

            prn_code(j,i) = prn1;           sv_num(j,i) = sv_num1;
            track_id(j,i) = find(trans_id_unique == transmitter_id1);

        end

    end
end

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

%% Part 3: L1a calibration
% this part converts from raw counts to signal power in watts and complete
% L1a calibration
clc

% use the upper most 5 delay rows as the noise floor
% TODO: still need to check the ddm_snr

% initialise variables for L1a results
power_analog = zeros(5,40,J,I)+invalid;
noise_floor = zeros(J,I)+invalid;
noise_power = zeros(J,I)+invalid;
snr_db = zeros(J,I)+invalid;

ddm_ant = zeros(J,I)+invalid;
inst_gain = zeros(J,I)+invalid;

for i = 1:I

    % retrieve noise standard deviation in counts for all three channels
    noise_std1 = [noise_std_dev_rf1(i),noise_std_dev_rf2(i),noise_std_dev_rf3(i)];

    for j = 1:J

        rf_source1 = rf_source(j,i);

        if ~isnan(rf_source1)

            first_scale_factor1 = first_scale_factor(j,i);
            ANZ_port1 = get_ANZ_port(rf_source1);
            raw_counts1 = raw_counts(:,:,j,i);                  % raw nadir antenna measurement in counts
            raw_counts1 = raw_counts1*first_scale_factor1;

            % perform L1a calibration
            %[signal_power_watts1,noise_power_watts1,noise_floor1] = L1a_counts2watts(raw_counts1, ...
            %    ANZ_port1,noise_std1);
            [signal_power_watts1,noise_power_watts1,noise_floor_counts1] = L1a_counts2watts1(raw_counts1,ANZ_port1);

            peak_raw_counts1 = max(raw_counts1,[],'all');
            peak_signal_watt1 = max(signal_power_watts1,[],'all');
            inst_gain1 = peak_raw_counts1/peak_signal_watt1;
            
            snr1 = peak_signal_watt1/noise_power_watts1;         
            snr1_db = pow2db(snr1);

            power_analog(:,:,j,i) = signal_power_watts1;
            noise_floor(j,i) = noise_floor_counts1;
            noise_power (j,i) = noise_power_watts1;
            snr_db(j,i) = snr1_db;

            inst_gain(j,i) = inst_gain1;
            ddm_ant(j,i) = ANZ_port1;            

        end
        
    end    
end

% save outputs to L1 structure
L1_postCal.raw_counts = raw_counts;
L1_postCal.l1a_power_ddm = power_analog;
L1_postCal.zenith_sig_i2q2 = zenith_i2q2;

L1_postCal.ddm_noise_floor = noise_floor;
L1_postCal.ddm_snr = snr_db;

L1_postCal.inst_gain = inst_gain;
L1_postCal.ddm_ant = ddm_ant;

% Part 3 ends

%% Part 4A: SP solver and geometries
clc

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

sx_rx_gain = zeros(J,I)+invalid;
%cross_pol = zeros(J,I)+invalid;

for i = 1:I

    % retrieve rx positions, velocities and attitdues
    rx_pos_xyz1 = rx_pos_xyz(i,:);      rx1.rx_pos_xyz = rx_pos_xyz1;
    rx_vel_xyz1 = rx_vel_xyz(i,:);      rx1.rx_vel_xyz = rx_vel_xyz1;
    rx_attitude1 = rx_attitude(i,:);    rx1.rx_attitude = rx_attitude1;    

    for j = 1:J

        % retrieve tx positions and velocities
        tx_pos_xyz1 = [tx_pos_x(j,i) tx_pos_y(j,i) tx_pos_z(j,i)];
        tx_vel_xyz1 = [tx_vel_x(j,i) tx_vel_y(j,i) tx_vel_z(j,i)];
        tx1.tx_pos_xyz = tx_pos_xyz1;
        tx1.tx_vel_xyz = tx_vel_xyz1;

        trans_id1 = prn_code(j,i);
        sv_num1 = sv_num(j,i);      tx1.sv_num = sv_num1;

        ddm_ant1 = ddm_ant(j,i);

        % only process these with valid svn code
        if sv_num1 ~= invalid

            % Part 4.1: SP solver
            % derive SP positions, angle of incidence and distance
            % to coast
            [sx_pos_xyz1,inc_angle_deg1,d_snell_deg1,dist_to_coast_km1,LOS_flag1] = sp_solver(tx_pos_xyz1,rx_pos_xyz1, ...
                dem,dtu10,landmask_nz);
            
            sx_pos_lla1 = ecef2lla(sx_pos_xyz1);            % <lat,lon,alt> of the specular reflection
            surface_type1 = get_surf_type(sx_pos_xyz1,landmask_nz,water_mask,lcv_mask);

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

            LOS_flag(j,i) = LOS_flag1;

            % define nadir antenna pattern to be used
            if ddm_ant1 == 2
                nadir_pattern = LHCP_pattern;

            elseif ddm_ant1 == 3
                nadir_pattern = RHCP_pattern;

            end

            % Part 4.2: SP-related variables - 1
            % this part derives tx/rx gains, ranges and other related
            % variables
            % only process samples with valid sx positions, i.e., LOS = 1
            if LOS_flag1 == 1
                
                % derive SP related geo-parameters, including angles
                % in various frames, ranges and antenna gain/GPS EIRP
                [sx_angle_body1,sx_angle_enu1,sx_angle_ant1,theta_gps1,ranges1,gps_rad1,rx_rad1] = spRelated(tx1,rx1, ...
                    sx_pos_xyz1,SV_eirp_LUT,nadir_pattern);

                % get values for deriving BRCS and reflectivity
                R_tsx1 = ranges1(1);        R_rsx1 = ranges1(2);
                gps_eirp_watt1 = gps_rad1(3);

                % derive cross polarisation
                if ddm_ant1 == 2
                    sx_rx_gain1 = rx_rad1(1);

                elseif ddm_ant1 == 3
                    sx_rx_gain1 = rx_rad1(2);

                end

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

                sx_rx_gain(j,i) = sx_rx_gain1;                  % nadir antenna gain - main polarisation
                              
            end
        end

    end
end

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

L1_postCal.sp_rx_gain = sx_rx_gain;

L1_postCal.gps_off_boresight_angle_deg = gps_boresight;

L1_postCal.static_gps_eirp = static_gps_eirp;
L1_postCal.gps_tx_power_db_w = gps_tx_power_db_w;
L1_postCal.gps_ant_gain_db_i = gps_ant_gain_db_i;

%% Part 4B: BRCS/NBRCS, reflectivity, coherent status and fresnel zone
clc

%load('../out/A_eff_large.mat')
%A_eff = A_eff_large;

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

norm_refl_waveform = zeros(1,40,J,I)+invalid;

nbrcs_scatter_area = zeros(J,I)+invalid;
ddm_nbrcs = zeros(J,I)+invalid;

les_scatter_area = zeros(J,I)+invalid;
ddm_les = zeros(J,I)+invalid;

tes_scatter_area = zeros(J,I)+invalid;
ddm_tes = zeros(J,I)+invalid;

surface_reflectivity = zeros(5,40,J,I)+invalid;
surface_reflectivity_peak = zeros(J,I)+invalid;

fresnel_coeff = zeros(J,I)+invalid;
fresnel_minor = zeros(J,I)+invalid;
fresnel_major = zeros(J,I)+invalid;
fresnel_orientation = zeros(J,I)+invalid;

coherency_ratio = zeros(J,I)+invalid;
coherency_state = zeros(J,I)+invalid;

for i = 1:I

    % retrieve rx positions and velocities
    rx_pos_xyz1 = rx_pos_xyz(i,:);          rx1.rx_pos_xyz = rx_pos_xyz1;
    rx_vel_xyz1 = rx_vel_xyz(i,:);          rx1.rx_vel_xyz = rx_vel_xyz1;
    rx_clk_drift1 = rx_clk_drift_mps(i,:);  rx1.rx_clk_drift = rx_clk_drift1;

    for j = 1:J

        % retrieve tx positions and velocities
        tx_pos_xyz1 = [tx_pos_x(j,i) tx_pos_y(j,i) tx_pos_z(j,i)];
        tx_vel_xyz1 = [tx_vel_x(j,i) tx_vel_y(j,i) tx_vel_z(j,i)];
        tx1.tx_pos_xyz = tx_pos_xyz1;
        tx1.tx_vel_xyz = tx_vel_xyz1;

        % retrieve sx-related parameters
        sx_pos_xyz1 = [sx_pos_x(j,i) sx_pos_y(j,i) sx_pos_z(j,i)];
        sx_pos_lla1 = ecef2lla(sx_pos_xyz1);

        sx_d_snell_deg1 = sx_d_snell_angle(j,i);
        dist_to_coast1 = dist_to_coast_km(j,i);

        sx1.sx_pos_xyz = sx_pos_xyz1;
        sx1.sx_d_snell = sx_d_snell_deg1;
        sx1.dist_to_coast = dist_to_coast1;

        eirp_watt1 = static_gps_eirp(j,i);
        rx_gain_db_i1 = sx_rx_gain(j,i);
        TSx1 = tx_to_sp_range(j,i);
        RSx1 = rx_to_sp_range(j,i);

        inc_angle1 = sx_inc_angle(j,i);
        ddm_ant1 = ddm_ant(j,i); 

        % retrieve ddm-related variables
        raw_counts1 = raw_counts(:,:,j,i);          ddm1.raw_counts = raw_counts1;
        power_analog1 = power_analog(:,:,j,i);      ddm1.power_analog1 = power_analog1;
        add_range_to_sp1 = add_range_to_sp(j,i);    ddm1.add_range_to_sp = add_range_to_sp1;
        snr_db1 = snr_db(j,i);                      ddm1.snr_db = snr_db1;
        
        delay_center_chips1 = delay_center_chips(j,i);
        ddm1.delay_center_chips = delay_center_chips1;

        doppler_center_hz1 = doppler_center_hz(j,i);
        ddm1.doppler_center_hz = doppler_center_hz1;

        T_coh1 = coherent_duration(i);
        ddm1.T_coh = T_coh1;

        ddm1.delay_resolution = 0.25;   ddm1.num_delay_bins = 40;   ddm1.delay_center_bin = 19;
        ddm1.doppler_resolution = 500;  ddm1.num_doppler_bins = 5;  ddm1.doppler_center_bin = 2;
        
        if (tx_pos_x(j,i) ~= invalid) && (sum(raw_counts1,'all')~=0)
            
            % Part 4.3: SP-related variables - 2
            % this part derives confidence, bin locations of SP
            raw_counts_max = max(raw_counts1,[],'all');
            [peak_doppler_bin1,peak_delay_bin1] = find(raw_counts1 == raw_counts_max,1);

            % TODO: correct doppler col and error to correct values once
            % center doppler values are correct
            [specular_bin1,zenith_code_phase1,confidence_flag1] = get_specular_bin(tx1,rx1,sx1,ddm1);
            
            % Part 4.4: brcs, nbrcs, effective area
            brcs1 = ddm_brcs(power_analog1,eirp_watt1,rx_gain_db_i1,TSx1,RSx1);

            %L = 9090; grid_res = 30; T_coh = 1/1000;
            %local_dem1 = get_local_dem(sx_pos_lla1,L,grid_res,dem,dtu10,dist_to_coast1);

            sx1.sx_delay_bin = specular_bin1(1);
            sx1.sx_doppler_bin = peak_doppler_bin1-1;             % TODO: change to SP dopp bin

            A_eff1 = ddm_Aeff(tx1,rx1,sx1,ddm1,local_dem1,T_coh);
            %A_eff1 = A_eff(:,:,j,i);
            [nbrcs1,LES1,TES1] = get_ddm_nbrcs(brcs1,A_eff1,sx1);
            
            % Part 4.5: reflectivity and peak reflectivity
            [refl1,refl_peak1] = ddm_refl(power_analog1,eirp_watt1,rx_gain_db_i1,TSx1,RSx1);

            % Part 4.6: Fresnel coefficient and dimensions
            [fresnel_coeff1,fresnel_axis1,fresnel_orientation1] = get_fresnel(tx1,rx1,sx1,inc_angle1,ddm_ant1);

            % Part 4.7: coherent status
            [CR1,CS1] = coh_det(raw_counts1,snr_db1);

            % normalised reflected waveform
            refl_waveform1 = sum(raw_counts1,1);
            norm_refl_waveform1 = refl_waveform1/max(refl_waveform1);
            
            % save to variables
            brcs_ddm_peak_bin_delay_row(j,i) = peak_delay_bin1-1;   % minus 1 for 0-based indces
            brcs_ddm_peak_bin_dopp_col(j,i) = peak_doppler_bin1-1;

            brcs_ddm_sp_bin_delay_row(j,i) = specular_bin1(1);
            brcs_ddm_sp_bin_dopp_col(j,i) = peak_doppler_bin1-1;    % TODO: specular_bin1(2);
            sp_delay_error(j,i) = specular_bin1(3);
            sp_dopp_error(j,i) = 0;                                 % TODO: specular_bin1(4);
            
            zenith_code_phase(j,i) = zenith_code_phase1;

            confidence_flag(j,i) = confidence_flag1;

            brcs(:,:,j,i) = brcs1;
            %A_eff(:,:,j,i) = A_eff1;
            
            nbrcs_scatter_area(j,i) = nbrcs1.nbrcs_scatter;
            ddm_nbrcs(j,i) = nbrcs1.nbrcs_value;

            les_scatter_area(j,i) = LES1.LES_scatter;
            ddm_les(j,i) = LES1.LES_slope;

            tes_scatter_area(j,i) = TES1.TES_scatter;
            ddm_tes(j,i) = TES1.TES_slope;

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

L1_postCal.brcs = brcs;
L1_postCal.eff_scatter = A_eff;

L1_postCal.nbrcs_scatter_area = nbrcs_scatter_area;
L1_postCal.ddm_nbrcs = ddm_nbrcs;

L1_postCal.les_scatter_area = les_scatter_area;
L1_postCal.ddm_les = ddm_les;

L1_postCal.tes_scatter_area = tes_scatter_area;
L1_postCal.ddm_tes = ddm_tes;

L1_postCal.surface_reflectivity = surface_reflectivity;
L1_postCal.surface_reflectivity_peak = surface_reflectivity_peak;

L1_postCal.fresnel_coeff = fresnel_coeff;
L1_postCal.fresnel_major = fresnel_major;
L1_postCal.fresnel_minor = fresnel_minor;
L1_postCal.fresnel_orientation = fresnel_orientation;

L1_postCal.coherency_ratio = coherency_ratio;
L1_postCal.coherency_state = coherency_state;

L1_postCal.norm_refl_waveform = norm_refl_waveform;

%% Cross Pol
clc

cross_pol = zeros(J,I)+invalid;

for i = 1:I
    for j = 1:J/2

        power_analog_L = power_analog(:,:,j,i);     max_power_L = max(max(power_analog_L));
        power_analog_R = power_analog(:,:,j+10,i);  max_power_R = max(max(power_analog_R));

        CP1 = max_power_L/max_power_R;
        CP_db1 = pow2db(CP1);

        cross_pol(j,i) = CP_db1;

    end
end

cross_pol(11:20,:) = -1*cross_pol(1:10,:);
L1_postCal.cross_pol = cross_pol;
L1_postCal.lna_noise_figure = zeros(J,I)+3;

%% Quality Flags
clc

quality_flags1 = zeros(J,I)+invalid;

for i = 1:I
    for j = 1:J

        quality_flag1_1 = zeros(1,23);

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
        signal_code_phase1 = delay_correction(meter2chips(add_range_to_sp(j,i)));
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
            det = find(brcs_ddma<0);
            if ~isempty(det)
                quality_flag1_1(17) = 1;
            end
        end

        % flag 18
        tx_pos_x1 = tx_pos_x(j,i);
        prn_code1 = prn_code(j,i);
        if (tx_pos_x1 == 0) && (prn_code1 ~= invalid)
            quality_flag1_1(18) = 1;
        end

        % flag 19
        sx_pos_x1 = sx_pos_x(j,i);
        if (sx_pos_x1 == invalid) && (prn_code1 ~= invalid)
            quality_flag1_1(19) = 1;
        end

        % flag 20
        rx_gain1 = sx_rx_gain(j,i);
        if (rx_gain1 == invalid) && (prn_code1 ~= invalid)
            quality_flag1_1(20) = 1;
        end

        quality_flag1_1(21) = 1;

        % flag 22
        rx_alt = rx_pos_lla(i,3);
        if rx_alt > 15000
            quality_flag1_1(22) = 1;
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
            quality_flag1_1(23) == 1 || ...
            quality_flag1_1(22) == 1 || ...
            quality_flag1_1(23) == 1)
        
            quality_flag1_1(1) = 1;

        end

        quality_flags1(j,i) = get_quality_flag(quality_flag1_1);
        
    end
end

L1_postCal.quality_flags1 = quality_flags1;

%% packet to netCDF
clc

netCDF_name = '../out/sample1.nc';
L1_dict = '../dat/L1_dict_final.xlsx';

L1_info = get_netcdf(netCDF_name,L1_dict,L1_postCal);

% L1 calibration ends

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function solves all L1 variables and packets as a structure for
% processing multiple L0 files
% Algorithm version 2.2
% L1 dictionary version 2.3

function L1_postCal = get_L1_product(   L0_filename, ...
                                        L1a_cal_ddm_counts_db,L1a_cal_ddm_power_dbm, ...
                                        dem,dtu10,landmask_nz,lcv_mask,water_mask, ...
                                        SV_PRN_LUT,SV_eirp_LUT,LHCP_pattern,RHCP_pattern, ...
                                        rx_alt_bins,inc_angle_bins,az_angle_bins,A_phy_LUT_all)

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

% delay and Doppler bin resolution
delay_bin_res = double(ncread(L0_filename,'/science/ddm/delay_bin_res_narrow'));
delay_bin_res = delay_bin_res(~isnan(delay_bin_res));
delay_bin_res = delay_bin_res(1);

doppler_bin_res = double(ncread(L0_filename,'/science/ddm/doppler_bin_res_narrow'));
doppler_bin_res = doppler_bin_res(~isnan(doppler_bin_res));
doppler_bin_res = doppler_bin_res(1);

% delay and Doppler center bin
center_delay_bin = double(ncread(L0_filename,'/science/ddm/ddm_center_delay_bin'));
center_delay_bin = center_delay_bin(~isnan(center_delay_bin));
center_delay_bin = center_delay_bin(1);

center_doppler_bin = double(ncread(L0_filename,'/science/ddm/ddm_center_doppler_bin'));
center_doppler_bin = center_doppler_bin(~isnan(center_doppler_bin));
center_doppler_bin = center_doppler_bin(1);

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

% DCP and FSW versions - new change 30 June
fsw_build_number = double(ncread(L0_filename,'/eng/fsw_build_number'));
fsw_major_version_number = double(ncread(L0_filename,'/eng/fsw_major_version_number'));
fsw_minor_version_number = double(ncread(L0_filename,'/eng/fsw_minor_version_number'));

fsw_version = [num2str(fsw_major_version_number(100)) '.' ...
               num2str(fsw_minor_version_number(100)) '.' ...
               num2str(fsw_build_number(100))];

dcp_build_number = double(ncread(L0_filename,'/eng/dcp_build_number'));
dcp_major_version_number = double(ncread(L0_filename,'/eng/dcp_major_version_number'));
dcp_minor_version_number = double(ncread(L0_filename,'/eng/dcp_minor_version_number'));

dcp_version = [num2str(dcp_major_version_number(100)) '.' ...
               num2str(dcp_minor_version_number(100)) '.' ...
               num2str(dcp_build_number(100))];

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

integration_duration = (coherent_duration.*non_coherent_integrations)*1000;

% temperatures from engineering data
index3 = ~isnan(eng_timestamp);

eng_timestamp = eng_timestamp(index3);
nadir_ant_temp_eng = nadir_ant_temp_eng(index3);
zenith_ant_temp_eng = zenith_ant_temp_eng(index3);

invalid = nan;                              % defines the value to be used for invalid fields

% M,N:number of doppler col and delay rows
% I,J:number of timestamps and NGRx capacity
[M,N,J,I] = size(raw_counts);

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
    int_duration1 = integration_duration(i);

    % convert pvt_gps_time to pvt_utc time
    D_pvt1 = gpstime2utc(pvt_gps_week1,pvt_gps_sec1);
    D_pvt2 = datetime(D_pvt1(1),D_pvt1(2),D_pvt1(3),D_pvt1(4),D_pvt1(5),D_pvt1(6)); % human time
    pvt_utc1 = convertTo(D_pvt2,'posixtime');    

    % derive ddm_utc, mid point of the non-coherent integrations
    ddm_utc1 = pvt_utc1+int_duration1/2;       
    
    % derive ddm_gps_week and ddm_gps_sec
    D_ddm1 = datetime(ddm_utc1,'ConvertFrom','posixtime');
    [year,month,date] = ymd(D_ddm1);
    [hour,mins,secs] = hms(D_ddm1);
    [ddm_gps_week1,ddm_gps_sec1] = utc2gpstime(year,month,date,hour,mins,secs);

    % save pvt and ddm timestamps for interpolation
    pvt_utc(i) = pvt_utc1;          ddm_utc(i) = ddm_utc1;
    gps_week(i) = ddm_gps_week1;    gps_tow(i) = ddm_gps_sec1;
    ddm_pvt_bias(i) = int_duration1/2;

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

% time coverage - updated compliance results 30 June
time_start = datetime(ddm_utc(1),'ConvertFrom','posixtime','format','yyyy-MM-dd HH:mm:ss');
time_start = char(time_start);
time_coverage_start = [time_start(1:10) 'T' time_start(end-7:end)];

time_end = datetime(ddm_utc(end),'ConvertFrom','posixtime','format','yyyy-MM-dd HH:mm:ss');
time_end = char(time_end);
time_coverage_end = [time_end(1:10) 'T' time_end(end-7:end)];

time_duration = ddm_utc(end)-ddm_utc(1)+1;
hours = floor(time_duration/3600);
minutes = floor((time_duration-hours*3600)/60);
seconds = time_duration-hours*3600-minutes*60;
time_coverage_duration = ['P0DT' num2str(hours) 'H' num2str(minutes) 'M' num2str(seconds) 'S'];

% write global variables - update compliance results 30 June
% reordering 30 June
L1_postCal.Conventions = 'CF-1.8, ACDD-1.3, ISO-8601';
L1_postCal.title = 'Rongowai Level 1 Science Data Record Version 1.0';
L1_postCal.history = join(['Mon Mar 20 22:21:12 2023: ncks -O -a -dsample,' ... % this line be the time produces L1
    '0,' num2str(I) ',1 -L1 --cnk_dmn=sample,1000' ....
    '--cnk_dmn=ddm,' num2str(J) '--cnk_dmn=delay,' num2str(M) '--cnk_dmn=doppler,' num2str(N) ...
    '--cnk_dmn=local_map_lat,41 --cnk_dmn=local_map_lon,41' ...
    '/tmp/qt_temp.MY8728 /tmp/qt_temp.jp8728\n./produce-L1-files' ...
    'production_1@rongowai-data-1.rongowai.auckland.ac.nz 7 ' ...
    time_coverage_start ' ' time_coverage_end ...
    ' --allow-partial --out-file=test-v1.nc']);
L1_postCal.standard_name_vocabulary = 'CF Standard Name Table v30';
L1_postCal.comment = 'DDMs are calibrated into Power (Watts) and Bistatic Radar Cross Section (m^2)';
L1_postCal.processing_level = '1';
L1_postCal.creator_type = 'institution';
L1_postCal.institution = 'University of Auckland (UoA)';
L1_postCal.creator_name = 'Rongowai Science Payloads Operations Centre';
L1_postCal.publisher_name = 'PO.DAAC';
L1_postCal.publisher_email = 'rongowai.auckland.ac.nz';
L1_postCal.publisher_url = 'spoc.auckland.ac.nz';
L1_postCal.geospatial_lat_min = '-48.034N';
L1_postCal.geospatial_lat_max = '-34.374N';
L1_postCal.geospatial_lon_min = '165.319E';
L1_postCal.geospatial_lon_max = '179.767E';

ref_timestamp_utc = ddm_utc(1);

pvt_timestamp_utc = pvt_utc-ref_timestamp_utc;
ddm_timestamp_utc = ddm_utc-ref_timestamp_utc;

L1_postCal.aircraft_reg = 'ZK-NFA';             % default value
L1_postCal.ddm_source = 2;                      % 1 = GPS signal simulator, 2 = aircraft
L1_postCal.ddm_time_type_selector = 1;          % 1 = middle of DDM sampling period
L1_postCal.delay_resolution = delay_bin_res;    % unit in chips
L1_postCal.dopp_resolution = doppler_bin_res;   % unit in Hz
L1_postCal.dem_source = 'SRTM30-200m';

% write algorithm and LUT versions
% version numbers may change when updating a LUT
L1_postCal.l1_algorithm_version = '2.2';        % 27 June, algorithm change
L1_postCal.l1_data_version = '2.3';             % 30 June, L1 dictionary change - compliance check
L1_postCal.l1a_sig_LUT_version = '1';
L1_postCal.l1a_noise_LUT_version = '1';
L1_postCal.A_LUT_version = '1.1';               % 27 June               
L1_postCal.ngrx_port_mapping_version = '1';
L1_postCal.nadir_ant_data_version = '1';
L1_postCal.zenith_ant_data_version = '1';
L1_postCal.nadir_ant_data_version = '2';        % 27 June
L1_postCal.prn_sv_maps_version = '1';
L1_postCal.gps_eirp_param_version = '7';
L1_postCal.land_mask_version = '1';
L1_postCal.surface_type_version = '1';
L1_postCal.mean_sea_surface_version = '1';
L1_postCal.per_bin_ant_version = '1';

L1_postCal.fsw_version = fsw_version;
L1_postCal.dcp_version = dcp_version;

L1_postCal.time_coverage_start = time_coverage_start;
L1_postCal.time_coverage_end = time_coverage_end;
L1_postCal.time_coverage_resolution = ddm_utc(2)-ddm_utc(1);
L1_postCal.time_coverage_duration = time_coverage_duration;

% write timestamps and ac-related variables
L1_postCal.pvt_timestamp_gps_week = pvt_gps_week;
L1_postCal.pvt_timestamp_gps_sec = pvt_gps_sec;
L1_postCal.pvt_timestamp_utc = pvt_timestamp_utc; 

L1_postCal.ddm_timestamp_gps_week = gps_week;
L1_postCal.ddm_timestamp_gps_sec = gps_tow;
L1_postCal.ddm_timestamp_utc = ddm_timestamp_utc;    
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
tx_pos_x(J/2+1:J,:) = tx_pos_x(1:J/2,:);    tx_vel_x(J/2+1:J,:) = tx_vel_x(1:J/2,:);        % XYZ error 27-June
tx_pos_y(J/2+1:J,:) = tx_pos_y(1:J/2,:);    tx_vel_y(J/2+1:J,:) = tx_vel_y(1:J/2,:);
tx_pos_z(J/2+1:J,:) = tx_pos_z(1:J/2,:);    tx_vel_z(J/2+1:J,:) = tx_vel_z(1:J/2,:);

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

% Part 3A: L1a calibration
% this part converts from raw counts to signal power in watts and complete
% L1a calibration

% initialise variables for L1a results
ddm_power_counts = zeros(M,N,J,I)+invalid;  % corrected raw counts with 1st sale factor
power_analog = zeros(M,N,J,I)+invalid;

ddm_ant = zeros(J,I)+invalid;
inst_gain = zeros(J,I)+invalid;

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
            % function update for this step 27-June
            ddm_power_watts1 = L1a_counts2watts(ddm_power_counts1,ANZ_port1, ...
                L1a_cal_ddm_counts_db,L1a_cal_ddm_power_dbm,std_dev1);

            % peak counts and power watts for instrument gain
            peak_counts1 = max(max(ddm_power_counts1));
            peak_power1 = max(max(ddm_power_watts1));

            inst_gain1 = peak_counts1/peak_power1;
            
            % save variables
            ddm_power_counts(:,:,j,i) = ddm_power_counts1;
            power_analog(:,:,j,i) = ddm_power_watts1;
            
            ddm_ant(j,i) = ANZ_port1;
            inst_gain(j,i) = inst_gain1;

        end
        
    end    
end

% save outputs to L1 structure
L1_postCal.raw_counts = ddm_power_counts;
L1_postCal.l1a_power_ddm = power_analog;
L1_postCal.zenith_sig_i2q2 = zenith_i2q2;

L1_postCal.inst_gain = inst_gain;
L1_postCal.ddm_ant = ddm_ant;

% Part 3A ends

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

% L1a confidence - 28 June
L1a_confidence_flag = zeros(J,I)+invalid;

for i = 1:I

    % retrieve rx positions, velocities and attitdues
    rx_pos_xyz1 = rx_pos_xyz(i,:);      rx1.rx_pos_xyz = rx_pos_xyz1;
    rx_vel_xyz1 = rx_vel_xyz(i,:);      rx1.rx_vel_xyz = rx_vel_xyz1;
    rx_attitude1 = [rx_attitude(i,1:2) 0];  % Euler angels are now in radians and yaw is resp. North   
    rx1.rx_attitude = rx_attitude1;    

    % variables are solved only for LHCP channels
    % RHCP channels share the same vales except RX gain solved for each
    % channel
    for j = 1:J/2

        % retrieve tx positions and velocities
        tx_pos_xyz1 = [tx_pos_x(j,i) tx_pos_y(j,i) tx_pos_z(j,i)];
        tx_vel_xyz1 = [tx_vel_x(j,i) tx_vel_y(j,i) tx_vel_z(j,i)];
        tx1.tx_pos_xyz = tx_pos_xyz1;
        tx1.tx_vel_xyz = tx_vel_xyz1;

        sv_num1 = sv_num(j,i);          tx1.sv_num = sv_num1;

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
                surface_type1 = get_surf_type2(sx_pos_xyz1,landmask_nz,water_mask,lcv_mask);

                % only process samples with valid sx positions, i.e., LOS = 1
                % derive sx velocity
                dt = 1;                                         % time step in second
                tx_pos_xyz_dt = tx_pos_xyz1+dt*tx_vel_xyz1;
                rx_pos_xyz_dt = rx_pos_xyz1+dt*rx_vel_xyz1;
                [sx_pos_xyz_dt,~,~,~,~] = sp_solver(tx_pos_xyz_dt,rx_pos_xyz_dt,dem,dtu10,landmask_nz);

                sx_vel_xyz1 = (sx_pos_xyz_dt-sx_pos_xyz1)/dt;            

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

                % determine L1a confidence - 28 June
                sx_theta_body1 = sx_angle_body1(1);         % off-boresight angle

                % antenna x-pol gain ratio
                copol_ratio1 = sx_rx_gain_LHCP1(1)-sx_rx_gain_LHCP1(2);
                xpol_ratio1 = sx_rx_gain_RHCP1(2)-sx_rx_gain_RHCP1(1);

                if sx_theta_body1<=60 && copol_ratio1>=14
                    L1a_confidence_flag_copol1 = 1;
                else
                    L1a_confidence_flag_copol1 = 0;
                end

                if sx_theta_body1<=60 && xpol_ratio1>=14
                    L1a_confidence_flag_xpol1 = 1;
                else
                    L1a_confidence_flag_xpol1 = 0;
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

                % copol gain
                sx_rx_gain_copol(j,i) = sx_rx_gain_LHCP1(1);      % LHCP channel LHCP gain
                sx_rx_gain_copol(j+J/2,i) = sx_rx_gain_RHCP1(2);  % RHCP channel RHCP gain

                % xpol gain
                sx_rx_gain_xpol(j,i) = sx_rx_gain_LHCP1(2);       % LHCP channel RHCP rx gain
                sx_rx_gain_xpol(j+J/2,i) = sx_rx_gain_RHCP1(1);   % RHCP channel LHCP rx gain

                % L1a confidence flag - 28 June
                L1a_confidence_flag(j,i) = L1a_confidence_flag_copol1;
                L1a_confidence_flag(j+J/2,i) = L1a_confidence_flag_xpol1;
                                              
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

% Part 4A ends

% Part 4B: Peak and theoretical SP index, zenith code phase 
% initialise variables
peak_delay_row = zeros(J,I)+invalid;
peak_doppler_col = zeros(J,I)+invalid;

sp_delay_row = zeros(J,I)+invalid;
sp_delay_error = zeros(J,I)+invalid;

sp_doppler_col = zeros(J,I)+invalid;
sp_doppler_error = zeros(J,I)+invalid;

zenith_code_phase = zeros(J,I)+invalid;

for i = 1:I

    rx_pos_xyz1 = rx_pos_xyz(i,:);
    rx_vel_xyz1 = rx_vel_xyz(i,:);

    for j = 1:J/2

        tx_pos_xyz1 = [tx_pos_x(j,i) tx_pos_y(j,i) tx_pos_z(j,i)];
        tx_vel_xyz1 = [tx_vel_x(j,i) tx_vel_y(j,i) tx_vel_z(j,i)];

        sx_pos_xyz1 = [sx_pos_x(j,i) sx_pos_y(j,i) sx_pos_z(j,i)];

        counts_LHCP1 = ddm_power_counts(:,:,j,i);
        
        add_range_to_sp1 = add_range_to_sp(j,i);        % from onboard tracker
        delay_center_chips1 = delay_center_chips(j,i);

        % zenith code phase
        add_range_to_sp_chips1 = meter2chips(add_range_to_sp1);
        zenith_code_phase1 = delay_center_chips1+add_range_to_sp_chips1;
        zenith_code_phase1 = delay_correction(zenith_code_phase1,1023);

        if ~isnan(tx_pos_x(j,i)) && ~isnan(sum(counts_LHCP1,'all'))

            % peak delay and doppler location
            % assume LHCP and RHCP DDMs have the same peak location
            peak_counts1 = max(max(counts_LHCP1));
            [peak_doppler_col1,peak_delay_row1] = find(counts_LHCP1==peak_counts1,1);

            % tx to rx range
            v_trx1 = tx_pos_xyz1-rx_pos_xyz1;   r_trx1 = norm(v_trx1);

            % SOC derived more accurate additional range to SP
            v_tsx1 = tx_pos_xyz1-sx_pos_xyz1;   r_tsx1 = norm(v_tsx1);
            v_rsx1 = rx_pos_xyz1-sx_pos_xyz1;   r_rsx1 = norm(v_rsx1);

            add_range_to_sp_soc1 = r_tsx1+r_rsx1-r_trx1;
            d_add_range1 = add_range_to_sp_soc1-add_range_to_sp1;

            d_delay_chips1 = meter2chips(d_add_range1);
            d_delay_bin1 = d_delay_chips1/delay_bin_res;
            
            sp_delay_row1 = center_delay_bin-d_delay_bin1;

            % SP doppler value
            [~,sp_doppler_hz1,~] = deldop(tx_pos_xyz1,rx_pos_xyz1, ...
                tx_vel_xyz1,rx_vel_xyz1,sx_pos_xyz1);

            doppler_center_hz1 = doppler_center_hz(j,i);

            % Doppler is now centrally binned - 27 June
            d_doppler_hz1 = doppler_center_hz1-sp_doppler_hz1;
            d_doppler_bin1 = d_doppler_hz1/doppler_bin_res;

            sp_doppler_col1 = center_doppler_bin-d_doppler_bin1;
            
            % SP delay and doppler location            
            peak_delay_row(j,i) = peak_delay_row1-1;            % correct to 0-based index
            peak_doppler_col(j,i) = peak_doppler_col1-1;

            sp_delay_row(j,i) = sp_delay_row1;
            sp_delay_error(j,i) = d_delay_chips1;

            sp_doppler_col(j,i) = sp_doppler_col1;
            sp_doppler_error(j,i) = d_doppler_hz1;

        end

        zenith_code_phase(j,i) = zenith_code_phase1;

    end
end

% extend to RHCP channels
peak_delay_row(J/2+1:J,:) = peak_delay_row(1:J/2,:);
peak_doppler_col(J/2+1:J,:) = peak_doppler_col(1:J/2,:);

sp_delay_row(J/2+1:J,:) = sp_delay_row(1:J/2,:);
sp_doppler_col(J/2+1:J,:) = sp_doppler_col(1:J/2,:);

% save variables
L1_postCal.brcs_ddm_peak_bin_delay_row = peak_delay_row;
L1_postCal.brcs_ddm_peak_bin_dopp_col = peak_doppler_col;

L1_postCal.brcs_ddm_sp_bin_delay_row = sp_delay_row;
L1_postCal.brcs_ddm_sp_bin_dopp_col = sp_doppler_col;

L1_postCal.sp_delay_error = sp_delay_error;
L1_postCal.sp_dopp_error = sp_doppler_error;

L1_postCal.zenith_code_phase = zenith_code_phase;

% Part 4B ends

% Part 3B: noise floor and SNR
delay_offset = 4;
sp_safe_margin = 9;         % safe space between SP and DDM end

noise_floor_all_LHCP = zeros(J/2,I)+invalid;
noise_floor_all_RHCP = zeros(J/2,I)+invalid;

% noise floor for all samples
for j = 1:J/2
    for i = I:I

        counts_LHCP1 = ddm_power_counts(:,:,j,i);
        counts_RHCP1 = ddm_power_counts(:,:,j+J/2,i);

        noise_floor_bins_LHCP1 = counts_LHCP1(:,end-delay_offset:end);
        noise_floor_bins_RHCP1 = counts_RHCP1(:,end-delay_offset:end);

        if ~isnan(tx_pos_x(j,i))

            noise_floor_LHCP1 = mean(noise_floor_bins_LHCP1,'all');
            noise_floor_all_LHCP(j,i) = noise_floor_LHCP1;

            noise_floor_RHCP1 = mean(noise_floor_bins_RHCP1,'all');
            noise_floor_all_RHCP(j,i) = noise_floor_RHCP1;

        end

    end
end

% single noise floor from valid DDMs
sp_delay_row_LHCP = sp_delay_row(1:10,:);       % reference to LHCP delay row

valid_idx = find(sp_delay_row_LHCP>0 & sp_delay_row_LHCP<39-sp_safe_margin & ...
    ~isnan(noise_floor_all_LHCP));

% noise floor is the median of the average counts
noise_floor_LHCP = median(noise_floor_all_LHCP(valid_idx));
noise_floor_RHCP = median(noise_floor_all_RHCP(valid_idx));

% SNR of SP
snr_LHCP_db = zeros(J/2,I)+invalid;
snr_flag_LHCP = zeros(J/2,I)+invalid;       % flag 0 for signal < 0

snr_RHCP_db = zeros(J/2,I)+invalid;
snr_flag_RHCP = zeros(J/2,I)+invalid;

for j = 1:J/2
    for i = 1:I

        counts_LHCP1 = ddm_power_counts(:,:,j,i);
        counts_RHCP1 = ddm_power_counts(:,:,j+J/2,i);

        sp_delay_row1 = floor(sp_delay_row_LHCP(j,i))+1;
        sp_doppler_col1 = floor(sp_doppler_col(j,i))+1;

        if sp_delay_row1<=40 && sp_delay_row1>0 && ...
                sp_doppler_col1<=5 && sp_doppler_col1>0

            sp_counts_LHCP1 = counts_LHCP1(sp_doppler_col1,sp_delay_row1);
            sp_counts_RHCP1 = counts_RHCP1(sp_doppler_col1,sp_delay_row1);

            signal_counts_LHCP1 = sp_counts_LHCP1-noise_floor_LHCP;
            snr_LHCP1 = signal_counts_LHCP1/noise_floor_LHCP;

            signal_counts_RHCP1 = sp_counts_RHCP1-noise_floor_RHCP;
            snr_RHCP1 = signal_counts_RHCP1/noise_floor_RHCP;

            if signal_counts_LHCP1>0
                snr_LHCP_db1 = pow2db(snr_LHCP1);
                snr_flag_LHCP1 = 1;

            elseif signal_counts_LHCP1<=0
                snr_LHCP_db1 = nan;
                snr_flag_LHCP1 = 0;

            end

            snr_LHCP_db(j,i) = snr_LHCP_db1;
            snr_flag_LHCP(j,i) = snr_flag_LHCP1;

            if signal_counts_RHCP1>0
                snr_RHCP_db1 = pow2db(snr_RHCP1);
                snr_flag_RHCP1 = 1;

            elseif signal_counts_RHCP1<=0
                snr_RHCP_db1 = nan;
                snr_flag_RHCP1 = 0;

            end

            snr_RHCP_db(j,i) = snr_RHCP_db1;
            snr_flag_RHCP(j,i) = snr_flag_RHCP1;            

        end

    end
end

% combine LHCP and RHCP results
noise_floor = [repmat(noise_floor_LHCP,[J/2,I]);repmat(noise_floor_RHCP,[J/2,I])];  % typo - 30 June
ddm_snr = [snr_LHCP_db;snr_RHCP_db];
snr_flag = [snr_flag_LHCP;snr_flag_RHCP];

L1_postCal.ddm_noise_floor = noise_floor;       % typo - 28 June
L1_postCal.ddm_snr = ddm_snr;
L1_postCal.ddm_snr_flag = snr_flag;             % typo - 28 June

% Part 3B ends

% Part 3C: confidence flag of the SP solved
confidence_flag = zeros(J,I)+invalid;

for j = 1:J/2
    for i = 1:I

        sx_delay_error1 = abs(sp_delay_error(j,i));
        sx_doppler_error1 = abs(sp_doppler_error(j,i));
        sx_d_snell_angle1 = abs(sx_d_snell_angle(j,i));

        snr_LHCP1 = snr_LHCP_db(j,i);

        if ~isnan(tx_pos_x(j,i))
    
            % criteria may change at a later stage
            delay_doppler_snell1 = (sx_delay_error1<1.25) && (abs(sx_doppler_error1)<250) && (sx_d_snell_angle1<2);
    
            if snr_LHCP1>=2 && ~delay_doppler_snell1
                confidence_flag1 = 0;
            elseif snr_LHCP1<2 && ~delay_doppler_snell1
                confidence_flag1 = 1;
            elseif snr_LHCP1<2 && delay_doppler_snell1
                confidence_flag1 = 2;
            elseif snr_LHCP1>=2 && delay_doppler_snell1
                confidence_flag1 = 3;
            else
                confidence_flag1 = nan;
            end

            confidence_flag(j,i) = confidence_flag1;

        end

    end
end

% expand to RHCP channels
confidence_flag(J/2+1:J,:) = confidence_flag(1:J/2,:);

L1_postCal.confidence_flag = confidence_flag;

L1_postCal.sp_ngrx_delay_correction = sp_delay_error;
L1_postCal.sp_ngrx_dopp_correction = sp_doppler_error;

% Part 5: Copol and xpol BRCS, reflectivity, peak reflectivity

% separate copol and xpol gain for using later
rx_gain_copol_LL = sx_rx_gain_copol(1:10,:);
rx_gain_copol_RR = sx_rx_gain_copol(11:20,:);

rx_gain_xpol_RL = sx_rx_gain_xpol(1:10,:);
rx_gain_xpol_LR = sx_rx_gain_xpol(11:20,:);

% BRCS, reflectivity 
brcs_copol = zeros(M,N,J/2,I)+invalid;
brcs_xpol  = zeros(M,N,J/2,I)+invalid;

refl_copol = zeros(M,N,J/2,I)+invalid;
refl_xpol  = zeros(M,N,J/2,I)+invalid;

sp_refl = zeros(J,I)+invalid;
norm_refl_waveform = zeros(1,40,J,I)+invalid;

for i = 1:I
    for j = 1:J/2

        cable_loss_db_LHCP = 0.6600;
        cable_loss_db_RHCP = 0.5840;

        % compensate cable loss
        power_analog_LHCP1 = power_analog(:,:,j,i)*db2pow(cable_loss_db_LHCP);
        power_analog_RHCP1 = power_analog(:,:,j+J/2,i)*db2pow(cable_loss_db_RHCP);

        R_tsx1 = tx_to_sp_range(j,i);
        R_rsx1 = rx_to_sp_range(j,i);

        rx_gain_dbi_1 = [   rx_gain_copol_LL(j,i) rx_gain_xpol_RL(j,i)  ;
                            rx_gain_xpol_LR(j,i)  rx_gain_copol_RR(j,i)];

        gps_eirp1 = static_gps_eirp(j,i);
        
        if ~isnan(sum(power_analog_LHCP1,'all'))
            [brcs_copol1,brcs_xpol1] = ddm_brcs2(power_analog_LHCP1,power_analog_RHCP1, ...
                gps_eirp1,rx_gain_dbi_1,R_tsx1,R_rsx1);

            [refl_copol1,refl_xpol1] = ddm_refl2(power_analog_LHCP1,power_analog_RHCP1, ...
                gps_eirp1,rx_gain_dbi_1,R_tsx1,R_rsx1);

            % reflectivity at SP
            sp_delay_row1 = floor(sp_delay_row(j,i))+1;
            sp_doppler_col1 = floor(sp_doppler_col(j,i))+1;

            if sp_delay_row1<=40 && sp_delay_row1>0 && ...
                    sp_doppler_col1<=5 && sp_doppler_col1>0
                sp_refl_copol1 = refl_copol1(sp_doppler_col1,sp_delay_row1);
                sp_refl_xpol1 = refl_xpol1(sp_doppler_col1,sp_delay_row1);

            else
                sp_refl_copol1 = nan;
                sp_refl_xpol1 = nan;

            end

            refl_waveform_copol1 = sum(refl_copol1,1);
            norm_refl_waveform_copol1 = refl_waveform_copol1/max(refl_waveform_copol1);

            refl_waveform_xpol1 = sum(refl_xpol1,1);
            norm_refl_waveform_xpol1 = refl_waveform_xpol1/max(refl_waveform_xpol1);

            brcs_copol(:,:,j,i) = brcs_copol1;
            brcs_xpol(:,:,j,i) = brcs_xpol1;

            refl_copol(:,:,j,i) = refl_copol1;
            refl_xpol(:,:,j,i) = refl_xpol1;

            sp_refl(j,i) = sp_refl_copol1;
            sp_refl(j+J/2,i) = sp_refl_xpol1;

            norm_refl_waveform(:,:,j,i) = norm_refl_waveform_copol1;
            norm_refl_waveform(:,:,j+J/2,i) = norm_refl_waveform_xpol1;

        end

    end
end

brcs(:,:,1:J/2,:) = brcs_copol;
brcs(:,:,J/2+1:J,:) = brcs_xpol;

surface_reflectivity(:,:,1:J/2,:) = refl_copol;
surface_reflectivity(:,:,J/2+1:J,:) = refl_copol;

L1_postCal.brcs = brcs;

L1_postCal.surface_reflectivity = surface_reflectivity;
L1_postCal.surface_reflectivity_peak = sp_refl;

L1_postCal.norm_refl_waveform = norm_refl_waveform;

% Part 5 ends

% Part 6: NBRCS and other related parameters
% significant changes in this section with a few new functions
% 28 June
% remove 3*3 NBRCS related scripts - 30 June
A_eff = zeros(M,N,J,I)+invalid;

nbrcs_scatter_area_v1 = zeros(J,I)+nan;
%nbrcs_scatter_area_v2 = zeros(J,I)+nan;

nbrcs_copol_v1 = zeros(J/2,I)+nan;
nbrcs_xpol_v1 = zeros(J/2,I)+nan;

%nbrcs_copol_v2 = zeros(J/2,I)+nan;
%nbrcs_xpol_v2 = zeros(J/2,I)+nan;

% 2D ambiguity function
chi2 = get_chi2(N,M,center_delay_bin+1,center_doppler_bin+1, ...
    delay_bin_res,doppler_bin_res);         % correct delay/Doppler index start from 1

for i = 1:I

    rx_vel_xyz1 = rx_vel_xyz(i,:);  rx1.rx_vel_xyz = rx_vel_xyz1;
    rx_alt1 = rx_pos_lla(i,3);

    for j = 1:J/2

        tx_vel_xyz1 = [tx_vel_x(j,i) tx_vel_y(j,i) tx_vel_z(j,i)];

        % azimuth angle between TX and RX velocity
        unit_rx_vel1 = rx_vel_xyz1/norm(rx_vel_xyz1);
        unit_tx_vel1 = tx_vel_xyz1/norm(tx_vel_xyz1);

        az_angle1 = acosd(-1*dot(unit_rx_vel1,unit_tx_vel1));   % 1st input of A_eff
        
        sx_pos_xyz1 = [sx_pos_x(j,i) sx_pos_y(j,i) sx_pos_z(j,i)];
        sx_lla1 = ecef2lla(sx_pos_xyz1);

        rx_alt_corrected1 = rx_alt1-sx_lla1(3);                 % 2nd input of A_eff

        inc_angle1 = sx_inc_angle(j,i);                         % 3rd input of A_eff
        
        brcs_copol1 = brcs_copol(:,:,j,i);
        brcs_xpol1 = brcs_xpol(:,:,j,i);

        %counts_LHCP1 = ddm_power_counts(:,:,j,i);
        %snr_LHCP1 = ddm_snr(j,i);
        
        % evaluate delay and Doppler bin location at SP
        sp_delay_row1 = sp_delay_row(j,i)+1;
        sp_doppler_col1 = sp_doppler_col(j,i)+1;

        % evaluate effective scattering area and NBRCS
        if sp_delay_row1<=39 && sp_delay_row1>=1 && ...
            sp_doppler_col1<=5 && sp_doppler_col1>=1 && ...     % secure the SP is within DDM range
            rx_alt_corrected1>=rx_alt_bins(1) && rx_alt_corrected1<=rx_alt_bins(end) && ...
            inc_angle1>=0 && inc_angle1<=80                     % ensure interpolate within reasonable range

            % effective scattering area
            A_eff1 = get_ddm_Aeff5(rx_alt1,inc_angle1,az_angle1, ...
                rx_alt_bins,inc_angle_bins,az_angle_bins, ...
                sp_delay_row1,sp_doppler_col1,chi2,A_phy_LUT_all);

            % nbrcs for SP bin
            [brcs_copol_ddma1,brcs_xpol_ddma1,A_eff_ddma1] = get_ddma_v1(brcs_copol1,brcs_xpol1,A_eff1, ...
                sp_delay_row1,sp_doppler_col1);

            nbrcs_copol_v1_1 = brcs_copol_ddma1/A_eff_ddma1;
            nbrcs_xpol_v1_1 = brcs_xpol_ddma1/A_eff_ddma1;

            nbrcs_scatter_area_v1(j,i) = A_eff_ddma1;
            
            nbrcs_copol_v1(j,i) = nbrcs_copol_v1_1;
            nbrcs_xpol_v1(j,i) = nbrcs_xpol_v1_1;

            % nbrcs for 3*3 bin
            %[brcs_copol_ddma2,brcs_xpol_ddma2,A_eff_ddma2] = get_ddma_v2(brcs_copol1,brcs_xpol1,A_eff1, ...
            %    sp_delay_row1,sp_doppler_col1);

            %nbrcs_copol_v2_1 = brcs_copol_ddma2/A_eff_ddma2;
            %nbrcs_xpol_v2_1 = brcs_xpol_ddma2/A_eff_ddma2;

            %nbrcs_scatter_area_v2(j,i) = A_eff_ddma2;
            
            %nbrcs_copol_v2(j,i) = nbrcs_copol_v2_1;
            %nbrcs_xpol_v2(j,i) = nbrcs_xpol_v2_1;
      
        end

    end
end

A_eff(:,:,J/2+1:J,:) = A_eff(:,:,1:J/2,:);

nbrcs_scatter_area_v1(J/2+1:J,:) = nbrcs_scatter_area_v1(1:J/2,:);
%nbrcs_scatter_area_v2(J/2+1:J,:) = nbrcs_scatter_area_v2(1:J/2,:);

ddm_nbrcs_v1 = [nbrcs_copol_v1;nbrcs_xpol_v1];
%ddm_nbrcs_v2 = [nbrcs_copol_v2;nbrcs_xpol_v2];

L1_postCal.eff_scatter = A_eff;

L1_postCal.nbrcs_scatter_area_v1 = nbrcs_scatter_area_v1;
%L1_postCal.nbrcs_scatter_area_v2 = nbrcs_scatter_area_v2;

L1_postCal.ddm_nbrcs_v1 = ddm_nbrcs_v1;
%L1_postCal.ddm_nbrcs_v2 = ddm_nbrcs_v2;

% Part 7: coherence detection
rmsd_delay_span_chips               = 1.5;  % delay span over which to compute error relative to WAF
fft_interpolation_factor            = 10;   % used in getRongowaiWAFRMSD for interpolation
power_vs_delay_noise_floor_rows     = 5;    % use in getRongowaiWAFRMSD power vs delay noise floor estimation
power_vs_delay_snr_min              = -10;  % minimum power vs delay snr threshold
aircraft_alt_min                    = 2e3;  % minimum height of aircraft, otherwise coherence state is 'uncertain'
max_peak_delay_row                  = 30;   % only use power vs delay waveforms where the peak is less than this

ac_alt = rx_pos_lla(:,3);

coherence_metric      = NaN(10,length(ac_alt));
coherence_state       = NaN(10,length(ac_alt));
power_vs_delay_snr_db = NaN(10,length(ac_alt));
peak_delay_row        = NaN(10,length(ac_alt));

for channel = 1:10      % channel-by-channel processing

    altitude  = ac_alt;
    ddm       = squeeze(ddm_power_counts(:,:,channel,:));     % select one channel at a time
    index     = 1:1:length(ac_alt);

    % select valid DDMs
    zeros_check = squeeze(sum(sum(ddm,1),2));           % select valid DDMs
    idx         = find(zeros_check>0);
    altitude    = altitude(idx);
    ddm         = ddm(:,:,idx);
    index       = index(idx);

    nan_check = squeeze(sum(sum(ddm,1),2));
    idx         = find(~isnan(nan_check));
    altitude    = altitude(idx);
    ddm         = ddm(:,:,idx);
    index       = index(idx);

    peak_delay_row_index = NaN(size(altitude));
    for z=1:length(altitude)
        i                           = index(z);
        hold_p_vs_delay             = squeeze(ddm(:,:,z));
        hold_p_vs_delay             = nansum(hold_p_vs_delay,1);
        [~,m]                       = nanmax(hold_p_vs_delay);
        peak_delay_row(channel,i)   = m;
        peak_delay_row_index(z)     = m;
    end

    idx             = find(peak_delay_row_index<=max_peak_delay_row); % detection not attempted if peak exceeds this limit
    ddm             = ddm(:,:,idx);
    index           = index(idx);
    altitude        = altitude(idx);

    for z=1:length(altitude)
        i = index(z);

        outputs                             = getRongowaiWAFRMSD(squeeze(ddm(:,:,z)),delay_bin_res, ...
                                                                 rmsd_delay_span_chips,fft_interpolation_factor, ...
                                                                 power_vs_delay_noise_floor_rows);

        coherence_metric(channel,i)         = outputs.rmsd;
        power_vs_delay_snr_db(channel,i)    = outputs.power_vs_delay_snr_db;

        if power_vs_delay_snr_db(channel,i)<power_vs_delay_snr_min || altitude(z)<aircraft_alt_min
            coherence_state(channel,i) = 5; % Uncertain coherence state
        elseif coherence_metric(channel,i)<=0.25
            coherence_state(channel,i) = 1; % With high confidence, state is dominantly coherent
        elseif coherence_metric(channel,i)>0.25 && coherence_metric(channel,i)<=0.50
            coherence_state(channel,i) = 2; % State is likely coherent
        elseif coherence_metric(channel,i)>0.50 && coherence_metric(channel,i)<0.75
            coherence_state(channel,i) = 3; % State is likely mixed/weakly diffuse
        elseif coherence_metric(channel,i)>=0.75
            coherence_state(channel,i) = 4; % With high confidence, state is dominantly incoherent
        end

    end 

end

coherence_metric(channel+1:channel+10,:) = coherence_metric;
coherence_state(channel+1:channel+10,:) = coherence_state;

L1_postCal.coherence_metric = coherence_metric;
L1_postCal.coherence_state = coherence_state;

% Part 8: fresnel dimensions
fresnel_coeff = zeros(J,I)+invalid;
fresnel_minor = zeros(J,I)+invalid;
fresnel_major = zeros(J,I)+invalid;
fresnel_orientation = zeros(J,I)+invalid;

for i = 1:I
    for j = 1:J

        tx_pos_xyz1 = [tx_pos_x(j,i) tx_pos_y(j,i) tx_pos_z(j,i)];
        rx_pos_xyz1 = rx_pos_xyz(i,:);
        sx_pos_xyz1 = [sx_pos_x(j,i) sx_pos_y(j,i) sx_pos_z(j,i)];

        inc_angle1 = sx_inc_angle(j,i);
        dist_to_coast1 = dist_to_coast_km(j,i);
        ddm_ant1 = ddm_ant(j,i);

        if ~isnan(ddm_ant1)

            [fresnel_coeff1,fresnel_axis1,fresnel_orientation1] = get_fresnel(tx_pos_xyz1, ...
                rx_pos_xyz1,sx_pos_xyz1,dist_to_coast1,inc_angle1,ddm_ant1);
    
            fresnel_coeff(j,i) = fresnel_coeff1;
            fresnel_major(j,i) = fresnel_axis1(1);
            fresnel_minor(j,i) = fresnel_axis1(2);
            fresnel_orientation(j,i) = fresnel_orientation1;

        end

    end
end

L1_postCal.fresnel_coeff = fresnel_coeff;
L1_postCal.fresnel_major = fresnel_major;
L1_postCal.fresnel_minor = fresnel_minor;
L1_postCal.fresnel_orientation = fresnel_orientation;

% Cross Pol - 28 June
% remove 3*3 related - 30 June
nbrcs_cross_pol_v1 = zeros(J,I)+invalid;
%nbrcs_cross_pol_v2 = zeros(J,I)+invalid;

for i = 1:I
    for j = 1:J/2

        nbrcs_LHCP1 = ddm_nbrcs_v1(j,i);
        nbrcs_RHCP1 = ddm_nbrcs_v1(j+J/2,i);

        %nbrcs_LHCP2 = ddm_nbrcs_v2(j,i);
        %nbrcs_RHCP2 = ddm_nbrcs_v2(j+J/2,i);

        CP1 = nbrcs_LHCP1/nbrcs_RHCP1;
        %CP2 = nbrcs_LHCP2/nbrcs_RHCP2;
        
        if CP1>0
            CP_db1 = pow2db(CP1);
            nbrcs_cross_pol_v1(j,i) = CP_db1;
            
        end

        %if CP2>0
        %    CP_db2 = pow2db(CP2);
        %    nbrcs_cross_pol_v2(j,i) = CP_db2;
            
        %end

    end
end

nbrcs_cross_pol_v1(11:20,:) = -1*nbrcs_cross_pol_v1(1:10,:);
%nbrcs_cross_pol_v2(11:20,:) = -1*nbrcs_cross_pol_v2(1:10,:);

L1_postCal.nbrcs_cross_pol_v1 = nbrcs_cross_pol_v1;
%L1_postCal.nbrcs_cross_pol_v2 = nbrcs_cross_pol_v2;

L1_postCal.lna_noise_figure = zeros(J,I)+3;         % LNA noise figure is 3 dB according to the specification

% Quality Flags
quality_flags1 = zeros(J,I)+invalid;

for i = 1:I
    for j = 1:J

        quality_flag1_1 = zeros(1,22);

        % flag 2, 3 and 4
        rx_roll1 = rad2deg(rx_roll(i));
        rx_pitch1 = rad2deg(rx_pitch(i));
        rx_yaw1 = rad2deg(rx_yaw(i));
    
        if abs(rx_roll1)<=30 && abs(rx_pitch1)<=10
            quality_flag1_1(2) = 1;
        else
            quality_flag1_1(3) = 1;
        end

        if rx_yaw1==-4
            quality_flag1_1(4) = 1;
        end

        % flag 5 - set if DDM is a test pattern - default 0
        quality_flag1_1(5) = 0;

        % flag 6 and 7
        trans_id1 = transmitter_id(j,i);
        if trans_id1 == 0
            quality_flag1_1(6) = 1;
        end

        if trans_id1 == 28
            quality_flag1_1(7) = 1;
        end

        % flag 8 and 11
        snr_db1 = ddm_snr(j,i);

        if i > 1
            snr_db2 = ddm_snr(j,i-1);
            diff1 = (db2pow(snr_db1)-db2pow(snr_db2))/db2pow(snr_db1);
            diff2 = snr_db1-snr_db2;

            if abs(diff1) > 0.1
                quality_flag1_1(7) = 1;
            end

            if abs(diff2) > 0.24
                quality_flag1_1(10) = 1;
            end
        end

        % flag 9 and 10
        dist_to_coast1 = dist_to_coast_km(j,i);

        if dist_to_coast1>0
            quality_flag1_1(9) = 1;
        end

        if dist_to_coast1>-5 && dist_to_caost1<0
            quality_flag1_1(10) = 1;
        end

        % flag 12
        ant_temp1 = ant_temp_nadir(i);
        if i > 1
            ant_temp2 = ant_temp_nadir(i-1);
            rate = (ant_temp2-ant_temp1)*60;

            if rate > 1
                quality_flag1_1(12) = 1;
            end

        end

        % flag 13
        zenith_code_phase1 = zenith_code_phase(j,i);
        signal_code_phase1 = delay_correction(meter2chips(add_range_to_sp(j,i)),1023);
        diff1 = zenith_code_phase1-signal_code_phase1;
        if diff1 >= 10
            quality_flag1_1(13) = 1;
        end

        % flag 14 and 15
        sp_delay_row1 = sp_delay_row(j,i);
        sp_dopp_col1 = sp_doppler_col(j,i);

        if (sp_delay_row1<14)||(sp_delay_row1>34)
            quality_flag1_1(14) = 1;
        end

        if (sp_dopp_col1<1)||(sp_dopp_col1>3)
            quality_flag1_1(15) = 1;
        end

        % flag 16
        ddma1 = nbrcs_scatter_area_v1(j,i);

        if isnan(ddma1) || ddma1<0
            quality_flag1_1(16) = 1;
        end

        % flag 17
        tx_pos_x1 = tx_pos_x(j,i);
        prn_code1 = prn_code(j,i);
        if (tx_pos_x1 == 0) && (~isnan(prn_code1))
            quality_flag1_1(17) = 1;
        end

        % flag 18
        sx_pos_x1 = sx_pos_x(j,i);
        if (sx_pos_x1 == invalid) && (~isnan(prn_code1))
            quality_flag1_1(18) = 1;
        end

        % flag 19
        rx_gain1 = sx_rx_gain_copol(j,i);
        if (rx_gain1 == invalid) && (~isnan(prn_code1))
            quality_flag1_1(19) = 1;
        end

        % flag 20 & 22
        rx_alt = rx_pos_lla(i,3);
        if rx_alt > 15000
            quality_flag1_1(20) = 1;
        end

        if rx_alt < 700
            quality_flag1_1(22) = 1;
        end

        % flag 21
        rx_vel_xyz1 = rx_vel_xyz(i,:);
        rx_speed1 = norm(rx_vel_xyz1);

        if rx_speed1>=150
            quality_flag1_1(21) = 1;
        end

        % flag 1
        if (quality_flag1_1(3) == 1 || ...
            quality_flag1_1(4) == 1 || ...
            quality_flag1_1(5) == 1 || ...
            quality_flag1_1(6) == 1 || ...
            quality_flag1_1(7) == 1 || ...
            quality_flag1_1(8) == 1 || ...
            quality_flag1_1(11) == 1 || ...
            quality_flag1_1(12) == 1 || ...
            quality_flag1_1(13) == 1 || ...
            quality_flag1_1(14) == 1 || ...
            quality_flag1_1(15) == 1 || ...
            quality_flag1_1(16) == 1 || ...
            quality_flag1_1(17) == 1 || ...
            quality_flag1_1(18) == 1 || ...
            quality_flag1_1(20) == 1 || ...
            quality_flag1_1(22) == 1)
        
            quality_flag1_1(1) = 1;

        end

        quality_flags1(j,i) = get_quality_flag(quality_flag1_1);
        
    end
end

%{
%=======
% This function solves all L1 variables and packets as a structure for
% processing multiple L0 files
% Algorithm version 2.2
% L1 dictionary version 2.2

function L1_postCal = get_L1_product(   L0_filename, ...
                                        L1a_cal_ddm_counts_db,L1a_cal_ddm_power_dbm, ...
                                        dem,dtu10,landmask_nz,lcv_mask,water_mask, ...
                                        SV_PRN_LUT,SV_eirp_LUT,LHCP_pattern,RHCP_pattern, ...
                                        rx_alt_bins,inc_angle_bins,az_angle_bins,A_phy_LUT_all)

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

% delay and Doppler bin resolution
delay_bin_res = double(ncread(L0_filename,'/science/ddm/delay_bin_res_narrow'));
delay_bin_res = delay_bin_res(~isnan(delay_bin_res));
delay_bin_res = delay_bin_res(1);

doppler_bin_res = double(ncread(L0_filename,'/science/ddm/doppler_bin_res_narrow'));
doppler_bin_res = doppler_bin_res(~isnan(doppler_bin_res));
doppler_bin_res = doppler_bin_res(1);

% delay and Doppler center bin
center_delay_bin = double(ncread(L0_filename,'/science/ddm/ddm_center_delay_bin'));
center_delay_bin = center_delay_bin(~isnan(center_delay_bin));
center_delay_bin = center_delay_bin(1);

center_doppler_bin = double(ncread(L0_filename,'/science/ddm/ddm_center_doppler_bin'));
center_doppler_bin = center_doppler_bin(~isnan(center_doppler_bin));
center_doppler_bin = center_doppler_bin(1);

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

% DCP and FSW versions - new change 30 June
fsw_build_number = double(ncread(L0_filename,'/eng/fsw_build_number'));
fsw_major_version_number = double(ncread(L0_filename,'/eng/fsw_major_version_number'));
fsw_minor_version_number = double(ncread(L0_filename,'/eng/fsw_minor_version_number'));

fsw_version = [num2str(fsw_major_version_number(100)) '.' ...
               num2str(fsw_minor_version_number(100)) '.' ...
               num2str(fsw_build_number(100))];

% DCP and FSW versions
dcp_build_number = double(ncread(L0_filename,'/eng/dcp_build_number'));
dcp_major_version_number = double(ncread(L0_filename,'/eng/dcp_major_version_number'));
dcp_minor_version_number = double(ncread(L0_filename,'/eng/dcp_minor_version_number'));

dcp_version = [num2str(dcp_major_version_number(100)) '.' ...
               num2str(dcp_minor_version_number(100)) '.' ...
               num2str(dcp_build_number(100))];

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

integration_duration = (coherent_duration.*non_coherent_integrations)*1000;

% temperatures from engineering data
index3 = ~isnan(eng_timestamp);

eng_timestamp = eng_timestamp(index3);
nadir_ant_temp_eng = nadir_ant_temp_eng(index3);
zenith_ant_temp_eng = zenith_ant_temp_eng(index3);

invalid = nan;                              % defines the value to be used for invalid fields

% M,N:number of doppler col and delay rows
% I,J:number of timestamps and NGRx capacity
[M,N,J,I] = size(raw_counts);

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
    int_duration1 = integration_duration(i);

    % convert pvt_gps_time to pvt_utc time
    D_pvt1 = gpstime2utc(pvt_gps_week1,pvt_gps_sec1);
    D_pvt2 = datetime(D_pvt1(1),D_pvt1(2),D_pvt1(3),D_pvt1(4),D_pvt1(5),D_pvt1(6)); % human time
    pvt_utc1 = convertTo(D_pvt2,'posixtime');    

    % derive ddm_utc, mid point of the non-coherent integrations
    ddm_utc1 = pvt_utc1+int_duration1/2;       
    
    % derive ddm_gps_week and ddm_gps_sec
    D_ddm1 = datetime(ddm_utc1,'ConvertFrom','posixtime');
    [year,month,date] = ymd(D_ddm1);
    [hour,mins,secs] = hms(D_ddm1);
    [ddm_gps_week1,ddm_gps_sec1] = utc2gpstime(year,month,date,hour,mins,secs);

    % save pvt and ddm timestamps for interpolation
    pvt_utc(i) = pvt_utc1;          ddm_utc(i) = ddm_utc1;
    gps_week(i) = ddm_gps_week1;    gps_tow(i) = ddm_gps_sec1;
    ddm_pvt_bias(i) = int_duration1/2;

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

% time coverage - updated compliance results 30 June
time_start = datetime(ddm_utc(1),'ConvertFrom','posixtime','format','yyyy-MM-dd HH:mm:ss');
time_start = char(time_start);
time_coverage_start = [time_start(1:10) 'T' time_start(end-7:end)];

time_end = datetime(ddm_utc(end),'ConvertFrom','posixtime','format','yyyy-MM-dd HH:mm:ss');
time_end = char(time_end);
time_coverage_end = [time_end(1:10) 'T' time_end(end-7:end)];

time_duration = ddm_utc(end)-ddm_utc(1)+1;
hours = floor(time_duration/3600);
minutes = floor((time_duration-hours*3600)/60);
seconds = time_duration-hours*3600-minutes*60;
time_coverage_duration = ['P0DT' num2str(hours) 'H' num2str(minutes) 'M' num2str(seconds) 'S'];

% write global variables - update compliance results 30 June
% reordering 30 June
L1_postCal.Conventions = 'CF-1.8, ACDD-1.3, ISO-8601';
L1_postCal.title = 'Rongowai Level 1 Science Data Record Version 1.0';
L1_postCal.history = join(['Mon Mar 20 22:21:12 2023: ncks -O -a -dsample,' ... % this line be the time produces L1
    '0,' num2str(I) ',1 -L1 --cnk_dmn=sample,1000' ....
    '--cnk_dmn=ddm,' num2str(J) '--cnk_dmn=delay,' num2str(M) '--cnk_dmn=doppler,' num2str(N) ...
    '--cnk_dmn=local_map_lat,41 --cnk_dmn=local_map_lon,41' ...
    '/tmp/qt_temp.MY8728 /tmp/qt_temp.jp8728\n./produce-L1-files' ...
    'production_1@rongowai-data-1.rongowai.auckland.ac.nz 7 ' ...
    time_coverage_start ' ' time_coverage_end ...
    ' --allow-partial --out-file=test-v1.nc']);
L1_postCal.standard_name_vocabulary = 'CF Standard Name Table v30';
L1_postCal.comment = 'DDMs are calibrated into Power (Watts) and Bistatic Radar Cross Section (m^2)';
L1_postCal.processing_level = '1';
L1_postCal.creator_type = 'institution';
L1_postCal.institution = 'University of Auckland (UoA)';
L1_postCal.creator_name = 'Rongowai Science Payloads Operations Centre';
L1_postCal.publisher_name = 'PO.DAAC';
L1_postCal.publisher_email = 'rongowai.auckland.ac.nz';
L1_postCal.publisher_url = 'spoc.auckland.ac.nz';
L1_postCal.geospatial_lat_min = '-48.034N';
L1_postCal.geospatial_lat_max = '-34.374N';
L1_postCal.geospatial_lon_min = '165.319E';
L1_postCal.geospatial_lon_max = '179.767E';

ref_timestamp_utc = ddm_utc(1);

pvt_timestamp_utc = pvt_utc-ref_timestamp_utc;
ddm_timestamp_utc = ddm_utc-ref_timestamp_utc;

L1_postCal.aircraft_reg = 'ZK-NFA';             % default value
L1_postCal.ddm_source = 2;                      % 1 = GPS signal simulator, 2 = aircraft
L1_postCal.ddm_time_type_selector = 1;          % 1 = middle of DDM sampling period
L1_postCal.delay_resolution = delay_bin_res;    % unit in chips
L1_postCal.dopp_resolution = doppler_bin_res;   % unit in Hz
L1_postCal.dem_source = 'SRTM30-200m';

% write algorithm and LUT versions
% version numbers may change when updating a LUT
L1_postCal.l1_algorithm_version = '2.2';        % 27 June, algorithm change
L1_postCal.l1_data_version = '2.3';             % 30 June, L1 dictionary change - compliance check
L1_postCal.l1a_sig_LUT_version = '1';
L1_postCal.l1a_noise_LUT_version = '1';
L1_postCal.A_LUT_version = '1.1';               % 27 June               
L1_postCal.ngrx_port_mapping_version = '1';
L1_postCal.nadir_ant_data_version = '1';
L1_postCal.zenith_ant_data_version = '1';
L1_postCal.nadir_ant_data_version = '2';        % 27 June
L1_postCal.prn_sv_maps_version = '1';
L1_postCal.gps_eirp_param_version = '7';
L1_postCal.land_mask_version = '1';
L1_postCal.surface_type_version = '1';
L1_postCal.mean_sea_surface_version = '1';
L1_postCal.per_bin_ant_version = '1';

L1_postCal.fsw_version = fsw_version;
L1_postCal.dcp_version = dcp_version;

L1_postCal.time_coverage_start = time_coverage_start;
L1_postCal.time_coverage_end = time_coverage_end;
L1_postCal.time_coverage_resolution = ddm_utc(2)-ddm_utc(1);
L1_postCal.time_coverage_duration = time_coverage_duration;

% write timestamps and ac-related variables
L1_postCal.pvt_timestamp_gps_week = pvt_gps_week;
L1_postCal.pvt_timestamp_gps_sec = pvt_gps_sec;
L1_postCal.pvt_timestamp_utc = pvt_timestamp_utc; 

L1_postCal.ddm_timestamp_gps_week = gps_week;
L1_postCal.ddm_timestamp_gps_sec = gps_tow;
L1_postCal.ddm_timestamp_utc = ddm_timestamp_utc;    
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
tx_pos_x(J/2+1:J,:) = tx_pos_x(1:J/2,:);    tx_vel_x(J/2+1:J,:) = tx_vel_x(1:J/2,:);        % XYZ error 27-June
tx_pos_y(J/2+1:J,:) = tx_pos_y(1:J/2,:);    tx_vel_y(J/2+1:J,:) = tx_vel_y(1:J/2,:);
tx_pos_z(J/2+1:J,:) = tx_pos_z(1:J/2,:);    tx_vel_z(J/2+1:J,:) = tx_vel_z(1:J/2,:);

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

% Part 3A: L1a calibration
% this part converts from raw counts to signal power in watts and complete
% L1a calibration

% initialise variables for L1a results
ddm_power_counts = zeros(M,N,J,I)+invalid;  % corrected raw counts with 1st sale factor
power_analog = zeros(M,N,J,I)+invalid;

ddm_ant = zeros(J,I)+invalid;
inst_gain = zeros(J,I)+invalid;

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
            % function update for this step 27-June
            ddm_power_watts1 = L1a_counts2watts(ddm_power_counts1,ANZ_port1, ...
                L1a_cal_ddm_counts_db,L1a_cal_ddm_power_dbm,std_dev1);

            % peak counts and power watts for instrument gain
            peak_counts1 = max(max(ddm_power_counts1));
            peak_power1 = max(max(ddm_power_watts1));

            inst_gain1 = peak_counts1/peak_power1;
            
            % save variables
            ddm_power_counts(:,:,j,i) = ddm_power_counts1;
            power_analog(:,:,j,i) = ddm_power_watts1;
            
            ddm_ant(j,i) = ANZ_port1;
            inst_gain(j,i) = inst_gain1;

        end
        
    end    
end

% save outputs to L1 structure
L1_postCal.raw_counts = ddm_power_counts;
L1_postCal.l1a_power_ddm = power_analog;
L1_postCal.zenith_sig_i2q2 = zenith_i2q2;

L1_postCal.inst_gain = inst_gain;
L1_postCal.ddm_ant = ddm_ant;

% Part 3A ends

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

% L1a confidence - 28 June
L1a_confidence_flag = zeros(J,I)+invalid;

for i = 1:I

    % retrieve rx positions, velocities and attitdues
    rx_pos_xyz1 = rx_pos_xyz(i,:);      rx1.rx_pos_xyz = rx_pos_xyz1;
    rx_vel_xyz1 = rx_vel_xyz(i,:);      rx1.rx_vel_xyz = rx_vel_xyz1;
    rx_attitude1 = [rx_attitude(i,1:2) 0];  % Euler angels are now in radians and yaw is resp. North   
    rx1.rx_attitude = rx_attitude1;    

    % variables are solved only for LHCP channels
    % RHCP channels share the same vales except RX gain solved for each
    % channel
    for j = 1:J/2

        % retrieve tx positions and velocities
        tx_pos_xyz1 = [tx_pos_x(j,i) tx_pos_y(j,i) tx_pos_z(j,i)];
        tx_vel_xyz1 = [tx_vel_x(j,i) tx_vel_y(j,i) tx_vel_z(j,i)];
        tx1.tx_pos_xyz = tx_pos_xyz1;
        tx1.tx_vel_xyz = tx_vel_xyz1;

        sv_num1 = sv_num(j,i);          tx1.sv_num = sv_num1;

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
                surface_type1 = get_surf_type2(sx_pos_xyz1,landmask_nz,water_mask,lcv_mask);

                % only process samples with valid sx positions, i.e., LOS = 1
                % derive sx velocity
                dt = 1;                                         % time step in second
                tx_pos_xyz_dt = tx_pos_xyz1+dt*tx_vel_xyz1;
                rx_pos_xyz_dt = rx_pos_xyz1+dt*rx_vel_xyz1;
                [sx_pos_xyz_dt,~,~,~,~] = sp_solver(tx_pos_xyz_dt,rx_pos_xyz_dt,dem,dtu10,landmask_nz);

                sx_vel_xyz1 = (sx_pos_xyz_dt-sx_pos_xyz1)/dt;            

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

                % determine L1a confidence - 28 June
                sx_theta_body1 = sx_angle_body1(1);         % off-boresight angle

                % antenna x-pol gain ratio
                copol_ratio1 = sx_rx_gain_LHCP1(1)-sx_rx_gain_LHCP1(2);
                xpol_ratio1 = sx_rx_gain_RHCP1(2)-sx_rx_gain_RHCP1(1);

                if sx_theta_body1<=60 && copol_ratio1>=14
                    L1a_confidence_flag_copol1 = 1;
                else
                    L1a_confidence_flag_copol1 = 0;
                end

                if sx_theta_body1<=60 && xpol_ratio1>=14
                    L1a_confidence_flag_xpol1 = 1;
                else
                    L1a_confidence_flag_xpol1 = 0;
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

                % copol gain
                sx_rx_gain_copol(j,i) = sx_rx_gain_LHCP1(1);      % LHCP channel LHCP gain
                sx_rx_gain_copol(j+J/2,i) = sx_rx_gain_RHCP1(2);  % RHCP channel RHCP gain

                % xpol gain
                sx_rx_gain_xpol(j,i) = sx_rx_gain_LHCP1(2);       % LHCP channel RHCP rx gain
                sx_rx_gain_xpol(j+J/2,i) = sx_rx_gain_RHCP1(1);   % RHCP channel LHCP rx gain

                % L1a confidence flag - 28 June
                L1a_confidence_flag(j,i) = L1a_confidence_flag_copol1;
                L1a_confidence_flag(j+J/2,i) = L1a_confidence_flag_xpol1;
                                              
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

% Part 4A ends

% Part 4B: Peak and theoretical SP index, zenith code phase 
% initialise variables
peak_delay_row = zeros(J,I)+invalid;
peak_doppler_col = zeros(J,I)+invalid;

sp_delay_row = zeros(J,I)+invalid;
sp_delay_error = zeros(J,I)+invalid;

sp_doppler_col = zeros(J,I)+invalid;
sp_doppler_error = zeros(J,I)+invalid;

zenith_code_phase = zeros(J,I)+invalid;

for i = 1:I

    rx_pos_xyz1 = rx_pos_xyz(i,:);
    rx_vel_xyz1 = rx_vel_xyz(i,:);

    for j = 1:J/2

        tx_pos_xyz1 = [tx_pos_x(j,i) tx_pos_y(j,i) tx_pos_z(j,i)];
        tx_vel_xyz1 = [tx_vel_x(j,i) tx_vel_y(j,i) tx_vel_z(j,i)];

        sx_pos_xyz1 = [sx_pos_x(j,i) sx_pos_y(j,i) sx_pos_z(j,i)];

        counts_LHCP1 = ddm_power_counts(:,:,j,i);
        
        add_range_to_sp1 = add_range_to_sp(j,i);        % from onboard tracker
        delay_center_chips1 = delay_center_chips(j,i);

        % zenith code phase
        add_range_to_sp_chips1 = meter2chips(add_range_to_sp1);
        zenith_code_phase1 = delay_center_chips1+add_range_to_sp_chips1;
        zenith_code_phase1 = delay_correction(zenith_code_phase1,1023);

        if ~isnan(tx_pos_x(j,i)) && ~isnan(sum(counts_LHCP1,'all'))

            % peak delay and doppler location
            % assume LHCP and RHCP DDMs have the same peak location
            peak_counts1 = max(max(counts_LHCP1));
            [peak_doppler_col1,peak_delay_row1] = find(counts_LHCP1==peak_counts1,1);

            % tx to rx range
            v_trx1 = tx_pos_xyz1-rx_pos_xyz1;   r_trx1 = norm(v_trx1);

            % SOC derived more accurate additional range to SP
            v_tsx1 = tx_pos_xyz1-sx_pos_xyz1;   r_tsx1 = norm(v_tsx1);
            v_rsx1 = rx_pos_xyz1-sx_pos_xyz1;   r_rsx1 = norm(v_rsx1);

            add_range_to_sp_soc1 = r_tsx1+r_rsx1-r_trx1;
            d_add_range1 = add_range_to_sp_soc1-add_range_to_sp1;

            d_delay_chips1 = meter2chips(d_add_range1);
            d_delay_bin1 = d_delay_chips1/delay_bin_res;
            
            sp_delay_row1 = center_delay_bin-d_delay_bin1;

            % SP doppler value
            [~,sp_doppler_hz1,~] = deldop(tx_pos_xyz1,rx_pos_xyz1, ...
                tx_vel_xyz1,rx_vel_xyz1,sx_pos_xyz1);

            doppler_center_hz1 = doppler_center_hz(j,i);

            % Doppler is now centrally binned - 27 June
            d_doppler_hz1 = doppler_center_hz1-sp_doppler_hz1;
            d_doppler_bin1 = d_doppler_hz1/doppler_bin_res;

            sp_doppler_col1 = center_doppler_bin-d_doppler_bin1;
            
            % SP delay and doppler location            
            peak_delay_row(j,i) = peak_delay_row1-1;            % correct to 0-based index
            peak_doppler_col(j,i) = peak_doppler_col1-1;

            sp_delay_row(j,i) = sp_delay_row1;
            sp_delay_error(j,i) = d_delay_chips1;

            sp_doppler_col(j,i) = sp_doppler_col1;
            sp_doppler_error(j,i) = d_doppler_hz1;

        end

        zenith_code_phase(j,i) = zenith_code_phase1;

    end
end

% extend to RHCP channels
peak_delay_row(J/2+1:J,:) = peak_delay_row(1:J/2,:);
peak_doppler_col(J/2+1:J,:) = peak_doppler_col(1:J/2,:);

sp_delay_row(J/2+1:J,:) = sp_delay_row(1:J/2,:);
sp_doppler_col(J/2+1:J,:) = sp_doppler_col(1:J/2,:);

% save variables
L1_postCal.brcs_ddm_peak_bin_delay_row = peak_delay_row;
L1_postCal.brcs_ddm_peak_bin_dopp_col = peak_doppler_col;

L1_postCal.brcs_ddm_sp_bin_delay_row = sp_delay_row;
L1_postCal.brcs_ddm_sp_bin_dopp_col = sp_doppler_col;

L1_postCal.sp_delay_error = sp_delay_error;
L1_postCal.sp_dopp_error = sp_doppler_error;

L1_postCal.zenith_code_phase = zenith_code_phase;

% Part 4B ends

% Part 3B: noise floor and SNR
delay_offset = 4;
sp_safe_margin = 9;         % safe space between SP and DDM end

noise_floor_all_LHCP = zeros(J/2,I)+invalid;
noise_floor_all_RHCP = zeros(J/2,I)+invalid;

% noise floor for all samples
for j = 1:J/2
    for i = I:I

        counts_LHCP1 = ddm_power_counts(:,:,j,i);
        counts_RHCP1 = ddm_power_counts(:,:,j+J/2,i);

        noise_floor_bins_LHCP1 = counts_LHCP1(:,end-delay_offset:end);
        noise_floor_bins_RHCP1 = counts_RHCP1(:,end-delay_offset:end);

        if ~isnan(tx_pos_x(j,i))

            noise_floor_LHCP1 = mean(noise_floor_bins_LHCP1,'all');
            noise_floor_all_LHCP(j,i) = noise_floor_LHCP1;

            noise_floor_RHCP1 = mean(noise_floor_bins_RHCP1,'all');
            noise_floor_all_RHCP(j,i) = noise_floor_RHCP1;

        end

    end
end

% single noise floor from valid DDMs
sp_delay_row_LHCP = sp_delay_row(1:10,:);       % reference to LHCP delay row

valid_idx = find(sp_delay_row_LHCP>0 & sp_delay_row_LHCP<39-sp_safe_margin & ...
    ~isnan(noise_floor_all_LHCP));

% noise floor is the median of the average counts
noise_floor_LHCP = median(noise_floor_all_LHCP(valid_idx));
noise_floor_RHCP = median(noise_floor_all_RHCP(valid_idx));

% SNR of SP
snr_LHCP_db = zeros(J/2,I)+invalid;
snr_flag_LHCP = zeros(J/2,I)+invalid;       % flag 0 for signal < 0

snr_RHCP_db = zeros(J/2,I)+invalid;
snr_flag_RHCP = zeros(J/2,I)+invalid;

for j = 1:J/2
    for i = 1:I

        counts_LHCP1 = ddm_power_counts(:,:,j,i);
        counts_RHCP1 = ddm_power_counts(:,:,j+J/2,i);

        sp_delay_row1 = floor(sp_delay_row_LHCP(j,i))+1;
        sp_doppler_col1 = floor(sp_doppler_col(j,i))+1;

        if sp_delay_row1<=40 && sp_delay_row1>0 && ...
                sp_doppler_col1<=5 && sp_doppler_col1>0

            sp_counts_LHCP1 = counts_LHCP1(sp_doppler_col1,sp_delay_row1);
            sp_counts_RHCP1 = counts_RHCP1(sp_doppler_col1,sp_delay_row1);

            signal_counts_LHCP1 = sp_counts_LHCP1-noise_floor_LHCP;
            snr_LHCP1 = signal_counts_LHCP1/noise_floor_LHCP;

            signal_counts_RHCP1 = sp_counts_RHCP1-noise_floor_RHCP;
            snr_RHCP1 = signal_counts_RHCP1/noise_floor_RHCP;

            if signal_counts_LHCP1>0
                snr_LHCP_db1 = pow2db(snr_LHCP1);
                snr_flag_LHCP1 = 1;

            elseif signal_counts_LHCP1<=0
                snr_LHCP_db1 = nan;
                snr_flag_LHCP1 = 0;

            end

            snr_LHCP_db(j,i) = snr_LHCP_db1;
            snr_flag_LHCP(j,i) = snr_flag_LHCP1;

            if signal_counts_RHCP1>0
                snr_RHCP_db1 = pow2db(snr_RHCP1);
                snr_flag_RHCP1 = 1;

            elseif signal_counts_RHCP1<=0
                snr_RHCP_db1 = nan;
                snr_flag_RHCP1 = 0;

            end

            snr_RHCP_db(j,i) = snr_RHCP_db1;
            snr_flag_RHCP(j,i) = snr_flag_RHCP1;            

        end

    end
end

% combine LHCP and RHCP results
noise_floor = [repmat(noise_floor_LHCP,[J/2,I]);repmat(noise_floor_RHCP,[J/2,I])];  % typo - 30 June
ddm_snr = [snr_LHCP_db;snr_RHCP_db];
snr_flag = [snr_flag_LHCP;snr_flag_RHCP];

L1_postCal.ddm_noise_floor = noise_floor;       % typo - 28 June
L1_postCal.ddm_snr = ddm_snr;
L1_postCal.ddm_snr_flag = snr_flag;             % typo - 28 June

% Part 3B ends

% Part 3C: confidence flag of the SP solved
confidence_flag = zeros(J,I)+invalid;

for j = 1:J/2
    for i = 1:I

        sx_delay_error1 = abs(sp_delay_error(j,i));
        sx_doppler_error1 = abs(sp_doppler_error(j,i));
        sx_d_snell_angle1 = abs(sx_d_snell_angle(j,i));

        snr_LHCP1 = snr_LHCP_db(j,i);

        if ~isnan(tx_pos_x(j,i))
    
            % criteria may change at a later stage
            delay_doppler_snell1 = (sx_delay_error1<1.25) && (abs(sx_doppler_error1)<250) && (sx_d_snell_angle1<2);
    
            if snr_LHCP1>=2 && ~delay_doppler_snell1
                confidence_flag1 = 0;
            elseif snr_LHCP1<2 && ~delay_doppler_snell1
                confidence_flag1 = 1;
            elseif snr_LHCP1<2 && delay_doppler_snell1
                confidence_flag1 = 2;
            elseif snr_LHCP1>=2 && delay_doppler_snell1
                confidence_flag1 = 3;
            else
                confidence_flag1 = nan;
            end

            confidence_flag(j,i) = confidence_flag1;

        end

    end
end

% expand to RHCP channels
confidence_flag(J/2+1:J,:) = confidence_flag(1:J/2,:);

L1_postCal.confidence_flag = confidence_flag;

L1_postCal.sp_ngrx_delay_correction = sp_delay_error;
L1_postCal.sp_ngrx_dopp_correction = sp_doppler_error;

% Part 5: Copol and xpol BRCS, reflectivity, peak reflectivity

% separate copol and xpol gain for using later
rx_gain_copol_LL = sx_rx_gain_copol(1:10,:);
rx_gain_copol_RR = sx_rx_gain_copol(11:20,:);

rx_gain_xpol_RL = sx_rx_gain_xpol(1:10,:);
rx_gain_xpol_LR = sx_rx_gain_xpol(11:20,:);

% BRCS, reflectivity 
brcs_copol = zeros(M,N,J/2,I)+invalid;
brcs_xpol  = zeros(M,N,J/2,I)+invalid;

refl_copol = zeros(M,N,J/2,I)+invalid;
refl_xpol  = zeros(M,N,J/2,I)+invalid;

sp_refl = zeros(J,I)+invalid;
norm_refl_waveform = zeros(1,40,J,I)+invalid;

for i = 1:I
    for j = 1:J/2

        cable_loss_db_LHCP = 0.6600;
        cable_loss_db_RHCP = 0.5840;

        % compensate cable loss
        power_analog_LHCP1 = power_analog(:,:,j,i)*db2pow(cable_loss_db_LHCP);
        power_analog_RHCP1 = power_analog(:,:,j+J/2,i)*db2pow(cable_loss_db_RHCP);

        R_tsx1 = tx_to_sp_range(j,i);
        R_rsx1 = rx_to_sp_range(j,i);

        rx_gain_dbi_1 = [   rx_gain_copol_LL(j,i) rx_gain_xpol_RL(j,i)  ;
                            rx_gain_xpol_LR(j,i)  rx_gain_copol_RR(j,i)];

        gps_eirp1 = static_gps_eirp(j,i);
        
        if ~isnan(sum(power_analog_LHCP1,'all'))
            [brcs_copol1,brcs_xpol1] = ddm_brcs2(power_analog_LHCP1,power_analog_RHCP1, ...
                gps_eirp1,rx_gain_dbi_1,R_tsx1,R_rsx1);

            [refl_copol1,refl_xpol1] = ddm_refl2(power_analog_LHCP1,power_analog_RHCP1, ...
                gps_eirp1,rx_gain_dbi_1,R_tsx1,R_rsx1);

            % reflectivity at SP
            sp_delay_row1 = floor(sp_delay_row(j,i))+1;
            sp_doppler_col1 = floor(sp_doppler_col(j,i))+1;

            if sp_delay_row1<=40 && sp_delay_row1>0 && ...
                    sp_doppler_col1<=5 && sp_doppler_col1>0
                sp_refl_copol1 = refl_copol1(sp_doppler_col1,sp_delay_row1);
                sp_refl_xpol1 = refl_xpol1(sp_doppler_col1,sp_delay_row1);

            else
                sp_refl_copol1 = nan;
                sp_refl_xpol1 = nan;

            end

            refl_waveform_copol1 = sum(refl_copol1,1);
            norm_refl_waveform_copol1 = refl_waveform_copol1/max(refl_waveform_copol1);

            refl_waveform_xpol1 = sum(refl_xpol1,1);
            norm_refl_waveform_xpol1 = refl_waveform_xpol1/max(refl_waveform_xpol1);

            brcs_copol(:,:,j,i) = brcs_copol1;
            brcs_xpol(:,:,j,i) = brcs_xpol1;

            refl_copol(:,:,j,i) = refl_copol1;
            refl_xpol(:,:,j,i) = refl_xpol1;

            sp_refl(j,i) = sp_refl_copol1;
            sp_refl(j+J/2,i) = sp_refl_xpol1;

            norm_refl_waveform(:,:,j,i) = norm_refl_waveform_copol1;
            norm_refl_waveform(:,:,j+J/2,i) = norm_refl_waveform_xpol1;

        end

    end
end

brcs(:,:,1:J/2,:) = brcs_copol;
brcs(:,:,J/2+1:J,:) = brcs_xpol;

surface_reflectivity(:,:,1:J/2,:) = refl_copol;
surface_reflectivity(:,:,J/2+1:J,:) = refl_copol;

L1_postCal.brcs = brcs;

L1_postCal.surface_reflectivity = surface_reflectivity;
L1_postCal.surface_reflectivity_peak = sp_refl;

L1_postCal.norm_refl_waveform = norm_refl_waveform;

% Part 5 ends

% Part 6: NBRCS and other related parameters
% significant changes in this section with a few new functions
% 28 June
A_eff = zeros(M,N,J,I)+invalid;

nbrcs_scatter_area_v1 = zeros(J,I)+nan;
nbrcs_scatter_area_v2 = zeros(J,I)+nan;

nbrcs_copol_v1 = zeros(J/2,I)+nan;
nbrcs_xpol_v1 = zeros(J/2,I)+nan;

nbrcs_copol_v2 = zeros(J/2,I)+nan;
nbrcs_xpol_v2 = zeros(J/2,I)+nan;

% 2D ambiguity function
chi2 = get_chi2(N,M,center_delay_bin+1,center_doppler_bin+1, ...
    delay_bin_res,doppler_bin_res);         % correct delay/Doppler index start from 1

for i = 1:I

    rx_vel_xyz1 = rx_vel_xyz(i,:);  rx1.rx_vel_xyz = rx_vel_xyz1;
    rx_alt1 = rx_pos_lla(i,3);

    for j = 1:J/2

        tx_vel_xyz1 = [tx_vel_x(j,i) tx_vel_y(j,i) tx_vel_z(j,i)];

        % azimuth angle between TX and RX velocity
        unit_rx_vel1 = rx_vel_xyz1/norm(rx_vel_xyz1);
        unit_tx_vel1 = tx_vel_xyz1/norm(tx_vel_xyz1);

        az_angle1 = acosd(-1*dot(unit_rx_vel1,unit_tx_vel1));   % 1st input of A_eff
        
        sx_pos_xyz1 = [sx_pos_x(j,i) sx_pos_y(j,i) sx_pos_z(j,i)];
        sx_lla1 = ecef2lla(sx_pos_xyz1);

        rx_alt_corrected1 = rx_alt1-sx_lla1(3);                 % 2nd input of A_eff

        inc_angle1 = sx_inc_angle(j,i);                         % 3rd input of A_eff
        
        brcs_copol1 = brcs_copol(:,:,j,i);
        brcs_xpol1 = brcs_xpol(:,:,j,i);

        counts_LHCP1 = ddm_power_counts(:,:,j,i);
        snr_LHCP1 = ddm_snr(j,i);
        
        % evaluate delay and Doppler bin location at SP
        sp_delay_row1 = sp_delay_row(j,i)+1;
        sp_doppler_col1 = sp_doppler_col(j,i)+1;

        % evaluate effective scattering area and NBRCS
        if sp_delay_row1<=39 && sp_delay_row1>=1 && ...
            sp_doppler_col1<=5 && sp_doppler_col1>=1 && ...     % secure the SP is within DDM range
            rx_alt_corrected1>=rx_alt_bins(1) && rx_alt_corrected1<=rx_alt_bins(end) && ...
            inc_angle1>=0 && inc_angle1<=80                     % ensure interpolate within reasonable range

            % effective scattering area
            A_eff1 = get_ddm_Aeff5(rx_alt1,inc_angle1,az_angle1, ...
                rx_alt_bins,inc_angle_bins,az_angle_bins, ...
                sp_delay_row1,sp_doppler_col1,chi2,A_phy_LUT_all);

            % nbrcs for SP bin
            [brcs_copol_ddma1,brcs_xpol_ddma1,A_eff_ddma1] = get_ddma_v1(brcs_copol1,brcs_xpol1,A_eff1, ...
                sp_delay_row1,sp_doppler_col1);

            nbrcs_copol_v1_1 = brcs_copol_ddma1/A_eff_ddma1;
            nbrcs_xpol_v1_1 = brcs_xpol_ddma1/A_eff_ddma1;

            nbrcs_scatter_area_v1(j,i) = A_eff_ddma1;
            
            nbrcs_copol_v1(j,i) = nbrcs_copol_v1_1;
            nbrcs_xpol_v1(j,i) = nbrcs_xpol_v1_1;

            % nbrcs for 3*3 bin
            [brcs_copol_ddma2,brcs_xpol_ddma2,A_eff_ddma2] = get_ddma_v1(brcs_copol1,brcs_xpol1,A_eff1, ...
                sp_delay_row1,sp_doppler_col1);

            nbrcs_copol_v2_1 = brcs_copol_ddma2/A_eff_ddma2;
            nbrcs_xpol_v2_1 = brcs_xpol_ddma2/A_eff_ddma2;

            nbrcs_scatter_area_v2(j,i) = A_eff_ddma2;
            
            nbrcs_copol_v2(j,i) = nbrcs_copol_v2_1;
            nbrcs_xpol_v2(j,i) = nbrcs_xpol_v2_1;
      
        end

    end
end

A_eff(:,:,J/2+1:J,:) = A_eff(:,:,1:J/2,:);

nbrcs_scatter_area_v1(J/2+1:J,:) = nbrcs_scatter_area_v1(1:J/2,:);
nbrcs_scatter_area_v2(J/2+1:J,:) = nbrcs_scatter_area_v2(1:J/2,:);

ddm_nbrcs_v1 = [nbrcs_copol_v1;nbrcs_xpol_v1];
ddm_nbrcs_v2 = [nbrcs_copol_v2;nbrcs_xpol_v2];

L1_postCal.eff_scatter = A_eff;

L1_postCal.nbrcs_scatter_area_v1 = nbrcs_scatter_area_v1;
L1_postCal.nbrcs_scatter_area_v2 = nbrcs_scatter_area_v2;

L1_postCal.ddm_nbrcs_v1 = ddm_nbrcs_v1;
L1_postCal.ddm_nbrcs_v2 = ddm_nbrcs_v2;

% Part 7: coherence detection
clc

rmsd_delay_span_chips               = 1.5;  % delay span over which to compute error relative to WAF
fft_interpolation_factor            = 10;   % used in getRongowaiWAFRMSD for interpolation
power_vs_delay_noise_floor_rows     = 5;    % use in getRongowaiWAFRMSD power vs delay noise floor estimation
power_vs_delay_snr_min              = -10;  % minimum power vs delay snr threshold
aircraft_alt_min                    = 2e3;  % minimum height of aircraft, otherwise coherence state is 'uncertain'
max_peak_delay_row                  = 30;   % only use power vs delay waveforms where the peak is less than this

ac_alt = rx_pos_lla(:,3);

coherence_metric      = NaN(10,length(ac_alt));
coherence_state       = NaN(10,length(ac_alt));
power_vs_delay_snr_db = NaN(10,length(ac_alt));
peak_delay_row        = NaN(10,length(ac_alt));

for channel = 1:10      % channel-by-channel processing

    altitude  = ac_alt;
    ddm       = squeeze(ddm_power_counts(:,:,channel,:));     % select one channel at a time
    index     = 1:1:length(ac_alt);

    % select valid DDMs
    zeros_check = squeeze(sum(sum(ddm,1),2));           % select valid DDMs
    idx         = find(zeros_check>0);
    altitude    = altitude(idx);
    ddm         = ddm(:,:,idx);
    index       = index(idx);

    nan_check = squeeze(sum(sum(ddm,1),2));
    idx         = find(~isnan(nan_check));
    altitude    = altitude(idx);
    ddm         = ddm(:,:,idx);
    index       = index(idx);

    peak_delay_row_index = NaN(size(altitude));
    for z=1:length(altitude)
        i                           = index(z);
        hold_p_vs_delay             = squeeze(ddm(:,:,z));
        hold_p_vs_delay             = nansum(hold_p_vs_delay,1);
        [~,m]                       = nanmax(hold_p_vs_delay);
        peak_delay_row(channel,i)   = m;
        peak_delay_row_index(z)     = m;
    end

    idx             = find(peak_delay_row_index<=max_peak_delay_row); % detection not attempted if peak exceeds this limit
    ddm             = ddm(:,:,idx);
    index           = index(idx);
    altitude        = altitude(idx);

    for z=1:length(altitude)
        i = index(z);

        outputs                             = getRongowaiWAFRMSD(squeeze(ddm(:,:,z)),delay_bin_res, ...
                                                                 rmsd_delay_span_chips,fft_interpolation_factor, ...
                                                                 power_vs_delay_noise_floor_rows);

        coherence_metric(channel,i)         = outputs.rmsd;
        power_vs_delay_snr_db(channel,i)    = outputs.power_vs_delay_snr_db;

        if power_vs_delay_snr_db(channel,i)<power_vs_delay_snr_min || altitude(z)<aircraft_alt_min
            coherence_state(channel,i) = 5; % Uncertain coherence state
        elseif coherence_metric(channel,i)<=0.25
            coherence_state(channel,i) = 1; % With high confidence, state is dominantly coherent
        elseif coherence_metric(channel,i)>0.25 && coherence_metric(channel,i)<=0.50
            coherence_state(channel,i) = 2; % State is likely coherent
        elseif coherence_metric(channel,i)>0.50 && coherence_metric(channel,i)<0.75
            coherence_state(channel,i) = 3; % State is likely mixed/weakly diffuse
        elseif coherence_metric(channel,i)>=0.75
            coherence_state(channel,i) = 4; % With high confidence, state is dominantly incoherent
        end

    end 

end

coherence_metric(channel+1:channel+10,:) = coherence_metric;
coherence_state(channel+1:channel+10,:) = coherence_state;

L1_postCal.coherence_metric = coherence_metric;
L1_postCal.coherence_state = coherence_state;

% Part 8: fresnel dimensions
fresnel_coeff = zeros(J,I)+invalid;
fresnel_minor = zeros(J,I)+invalid;
fresnel_major = zeros(J,I)+invalid;
fresnel_orientation = zeros(J,I)+invalid;

for i = 1:I
    for j = 1:J

        tx_pos_xyz1 = [tx_pos_x(j,i) tx_pos_y(j,i) tx_pos_z(j,i)];
        rx_pos_xyz1 = rx_pos_xyz(i,:);
        sx_pos_xyz1 = [sx_pos_x(j,i) sx_pos_y(j,i) sx_pos_z(j,i)];

        inc_angle1 = sx_inc_angle(j,i);
        dist_to_coast1 = dist_to_coast_km(j,i);
        ddm_ant1 = ddm_ant(j,i);

        if ~isnan(ddm_ant1)

            [fresnel_coeff1,fresnel_axis1,fresnel_orientation1] = get_fresnel(tx_pos_xyz1, ...
                rx_pos_xyz1,sx_pos_xyz1,dist_to_coast1,inc_angle1,ddm_ant1);
    
            fresnel_coeff(j,i) = fresnel_coeff1;
            fresnel_major(j,i) = fresnel_axis1(1);
            fresnel_minor(j,i) = fresnel_axis1(2);
            fresnel_orientation(j,i) = fresnel_orientation1;

        end

    end
end

L1_postCal.fresnel_coeff = fresnel_coeff;
L1_postCal.fresnel_major = fresnel_major;
L1_postCal.fresnel_minor = fresnel_minor;
L1_postCal.fresnel_orientation = fresnel_orientation;

% Cross Pol - 28 June
nbrcs_cross_pol_v1 = zeros(J,I)+invalid;
nbrcs_cross_pol_v2 = zeros(J,I)+invalid;

for i = 1:I
    for j = 1:J/2

        nbrcs_LHCP1 = ddm_nbrcs_v1(j,i);
        nbrcs_RHCP1 = ddm_nbrcs_v1(j+J/2,i);

        nbrcs_LHCP2 = ddm_nbrcs_v2(j,i);
        nbrcs_RHCP2 = ddm_nbrcs_v2(j+J/2,i);

        CP1 = nbrcs_LHCP1/nbrcs_RHCP1;
        CP2 = nbrcs_LHCP2/nbrcs_RHCP2;
        
        if CP1>0
            CP_db1 = pow2db(CP1);
            nbrcs_cross_pol_v1(j,i) = CP_db1;
            
        end

        if CP2>0
            CP_db2 = pow2db(CP2);
            nbrcs_cross_pol_v2(j,i) = CP_db2;
            
        end

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

        quality_flag1_1 = zeros(1,26);

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
        snr_db1 = ddm_snr(j,i);

        if i > 1
            snr_db2 = ddm_snr(j,i-1);
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
        sp_delay_row1 = sp_delay_row(j,i);
        sp_dopp_col1 = sp_doppler_col(j,i);

        if (sp_delay_row1<15)||(sp_delay_row1>35)
            quality_flag1_1(15) = 1;
        end

        if (sp_dopp_col1<2)||(sp_dopp_col1>4)
            quality_flag1_1(16) = 1;
        end

        % flag 17
        if (floor(sp_delay_row1) < 38) && (floor(sp_delay_row1) > 0) && ...
                (floor(sp_dopp_col1) < 5) && (floor(sp_dopp_col1) > 1)
            brcs_ddma = brcs(floor(sp_dopp_col1)-1:floor(sp_dopp_col1)+1, ...
                            floor(sp_delay_row1):floor(sp_dopp_col1)+3);
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

        % flag 22 & 26
        rx_alt = rx_pos_lla(i,3);
        if rx_alt > 15000
            quality_flag1_1(22) = 1;
        end

        if rx_alt < 700
            quality_flag1_1(26) = 1;
        end

        % flag 24
        prn1 = prn_code(j,i);
        if prn1 == 28
            quality_flag1_1(24) = 1;
        end

        % flag 25
        rx_vel_xyz1 = rx_vel_xyz(i,:);
        rx_speed1 = norm(rx_vel_xyz1);

        if rx_speed1>150
            quality_flag1_1(25) = 1;
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
            quality_flag1_1(24) == 1 || ...
            quality_flag1_1(25) == 1 || ...
            quality_flag1_1(26) == 1)
        
            quality_flag1_1(1) = 1;

        end

        quality_flags1(j,i) = get_quality_flag(quality_flag1_1);
        
    end
end

>>>>>>> 533a4a80b6c994674a6e1315cb442a711c900a39
L1_postCal.quality_flags1 = quality_flags1;
%}
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
%L0_filename = [path1 '2022-09-22_L0.nc'];
L0_filename = [path1 'Flight_1_2022-10-02.nc'];

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
rx_pitch_deg_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_attitude_pitch_deg'));
rx_roll_deg_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_attitude_roll_deg'));
rx_yaw_deg_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_attitude_yaw_deg'));

% rx clock bias and drifts
rx_clk_bias_m_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_clock_bias_m'));
rx_clk_drift_mps_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_clock_drift_mps'));

% tx ID/satellite PRN
transmitter_id = double(ncread(L0_filename,'/science/ddm/transmitter_id'));

% raw counts and ddm parameters
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
doppler_center_Hz = double(ncread(L0_filename,'/science/ddm/center_doppler_bin_frequency'));

% number of doppler and delay bins
num_delay_bins = double(ncread(L0_filename,'/science/ddm/num_delay_bins')); 
num_doppler_bins = double(ncread(L0_filename,'/science/ddm/num_doppler_bins'));

% coherent duration and noncoherent integration
% TODO: double check the units
coherent_duration = double(ncread(L0_filename,'/science/ddm/L1_E1_coherent_duration'));
non_coherent_integrations = double(ncread(L0_filename,'/science/ddm/L1_E1_non_coherent_integrations'));

coherent_duration = coherent_duration/1000;
non_coherent_integrations = non_coherent_integrations/1000;

% NGRx estimate additional delay path
add_range_to_sp = double(ncread(L0_filename,'/science/ddm/additional_range_to_SP'));

wgs84_flag = 0;             % only for bench test data

% Load L0 data ends

%% Prelaunch 2 - define external data paths and filenames
clc

% define IGS orbits filename (*.sp3)
% note this path is defined in C++ format
gps_orbit_filename = '..//dat//orbits//igr21526.sp3';

% load SRTM_30 DEM
dem_path = '../dat/dem/';
dem_file1 = 'nzsrtm_30_part1_v1.dat';
dem_file2 = 'nzsrtm_30_part2_v1.dat';

dem1 = get_dem([dem_path dem_file1]);
dem2 = get_dem([dem_path dem_file2]);

if wgs84_flag == 1
    dem1.ele = zeros(length(dem1.lat),length(dem1.lon));
    dem2.ele = zeros(length(dem2.lat),length(dem2.lon));
end

% load DTU10 model
dtu_path = '../dat/dtu/';
dtu_filename = 'dtu10_v1.dat';

dtu10 = get_dtu10([dtu_path dtu_filename]);

if wgs84_flag == 1
    dtu10.ele = zeros(length(dtu10.lat),length(dtu10.lon));
end

% load ocean/land (distance to coast) mask
landmask_path = '../dat/cst/';
landmask_filename = 'dist_to_coast_nz_v1.dat';

landmask_nz = get_dist_to_coast_mask([landmask_path landmask_filename]);

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

%% Prelaunch 2.5: Filter valid timestampes
clc

valid_index = ~isnan(pvt_gps_week);

pvt_gps_week = pvt_gps_week(valid_index);       pvt_gps_sec = pvt_gps_sec(valid_index);

rx_pos_x_pvt = rx_pos_x_pvt(valid_index);       rx_pos_y_pvt = rx_pos_y_pvt(valid_index);
rx_pos_z_pvt = rx_pos_z_pvt(valid_index);

rx_vel_x_pvt = rx_vel_x_pvt(valid_index);       rx_vel_y_pvt = rx_vel_y_pvt(valid_index);
rx_vel_z_pvt = rx_vel_z_pvt(valid_index);

rx_roll_deg_pvt = rx_roll_deg_pvt(valid_index); rx_pitch_deg_pvt = rx_pitch_deg_pvt(valid_index);
rx_yaw_deg_pvt = rx_yaw_deg_pvt(valid_index);

rx_clk_bias_m_pvt = rx_clk_bias_m_pvt(valid_index);
rx_clk_drift_mps_pvt = rx_clk_drift_mps_pvt(valid_index);

transmitter_id = transmitter_id(:,valid_index);

raw_counts = raw_counts(:,:,:,valid_index);     
zenith_i2q2 = zenith_i2q2(:,valid_index);

rf_source = rf_source(:,valid_index);
ddm_number = ddm_number(:,valid_index);

noise_std_dev_rf1 = noise_std_dev_rf1(valid_index);
noise_std_dev_rf2 = noise_std_dev_rf2(valid_index); 
noise_std_dev_rf3 = noise_std_dev_rf3(valid_index);

delay_bin_res = delay_bin_res(valid_index);     doppler_bin_res = doppler_bin_res(valid_index);

delay_center_bin = delay_center_bin(valid_index);
doppler_center_bin = doppler_center_bin(valid_index);

% absolute ddm center delay and doppler
delay_center_chips = delay_center_chips(:,valid_index);    
doppler_center_Hz = doppler_center_Hz(:,valid_index);

% number of doppler and delay bins
num_delay_bins = num_delay_bins(valid_index);   num_doppler_bins = num_doppler_bins(valid_index);

% coherent duration and noncoherent integration
% TODO: double check the units
coherent_duration = coherent_duration(valid_index);
non_coherent_integrations = non_coherent_integrations(valid_index);

%coherent_duration = coherent_duration/1000;
%non_coherent_integrations = non_coherent_integrations/1000;

% NGRx estimate additional delay path
add_range_to_sp = add_range_to_sp(:,valid_index);

%% Part 1: General processing
% 1 - process timestamps and UTC human time
% 2 - interpolate all parameters to DDM timestamp
% 3 - other general/global variables
clc

invalid = -9999;                                    % defines the value to be used for invalid fields

I = length(pvt_gps_sec);                            % total length of samples
J = 20;                                             % maximal NGRx capacity

non_coherent_integrations(isnan(non_coherent_integrations)) = invalid;

% initialise output data array for Part 1 processing
pvt_utc = zeros(I,1)+invalid;       ddm_utc = zeros(I,1)+invalid;
gps_week = zeros(I,1)+invalid;      gps_tow = zeros(I,1)+invalid;
ddm_pvt_bias = zeros(I,1)+invalid;

% initialise a structure to save L1 results
L1_postCal = struct;

% derive and save ddm_timestamp_utc/gps to L1 structure
for i = 1:I

    pvt_gps_week1 = pvt_gps_week(i);
    pvt_gps_sec1 = pvt_gps_sec(i);
    non_coherent1 = non_coherent_integrations(i);

    if non_coherent1 ~= -9999
     
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

end

% linear interpolation all the values at ddm timestamp
[rx_pos_xyz,rx_vel_xyz,rx_attitude,rx_clk] = PVT2ddm_timestamp(pvt_utc,ddm_utc,rx_pos_xyz_pvt, ...
    rx_vel_xyz_pvt,rx_attitude_pvt,rx_clk_pvt);

% get <lat,lon,alt> of aircraft
rx_pos_lla = ecef2lla(rx_pos_xyz);

% write global variables
% TODO: flight ID needs to be retrieved from OpenSky
L1_postCal.time_coverage_start = datetime(ddm_utc(1),'ConvertFrom','posixtime');        % human time
L1_postCal.time_coverage_end = datetime(ddm_utc(end),'ConvertFrom','posixtime');
L1_postCal.time_coverage_resolution = ddm_utc(2)-ddm_utc(1);

% time coverage, not possible for more than 1 day
time_duration = ddm_utc(end)-ddm_utc(1)+1;
hours = floor(time_duration/3600);
minutes = floor((time_duration-hours*3600)/60);
seconds = time_duration-hours*3600-minutes*60;
time_coverage_duration = ['P0DT' num2str(hours) 'H' num2str(minutes) 'M' num2str(seconds) 'S'];

L1_postCal.time_coverage_duration = time_coverage_duration;

L1_postCal.aircraft_reg = 'ZK-NFA';
L1_postCal.ddm_source = 1;                      % 1 = GPS signal simulator, 2 = aircraft
L1_postCal.ddm_time_type_selector = 1;          % 1 = middle of DDM sampling period
L1_postCal.delay_resolution = 0.25;             % unit in chips
L1_postCal.dopp_resolution = 500;               % unit in Hz
L1_postCal.fixed_noise_power = 78.3;            % unit in dB
L1_postCal.dem_source = 'SRTM30';

% write version numbers
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
L1_postCal.mean_sea_surface_version = '1';

% write timestamps and ac-related variables
L1_postCal.pvt_timestamp_gps_week = pvt_gps_week;
L1_postCal.pvt_timestamp_gps_sec = pvt_gps_sec;
L1_postCal.pvt_timestamp_utc = pvt_utc;

L1_postCal.ddm_timestamp_gps_week = gps_week;
L1_postCal.ddm_timestamp_gps_sec = gps_tow;
L1_postCal.ddm_timestamp_utc = ddm_utc;

L1_postCal.ddm_pvt_bias = ddm_pvt_bias;

L1_postCal.ac_lat = rx_pos_lla(:,1);
L1_postCal.ac_lon = rx_pos_lla(:,2);
L1_postCal.ac_alt = rx_pos_lla(:,3);

L1_postCal.ac_pos_x_pvt = rx_pos_x_pvt;
L1_postCal.ac_pos_y_pvt = rx_pos_y_pvt;
L1_postCal.ac_pos_z_pvt = rx_pos_z_pvt;

L1_postCal.ac_pos_x = rx_pos_xyz(:,1);
L1_postCal.ac_pos_y = rx_pos_xyz(:,2);
L1_postCal.ac_pos_z = rx_pos_xyz(:,3);

L1_postCal.ac_vel_x_pvt = rx_vel_x_pvt;
L1_postCal.ac_vel_y_pvt = rx_vel_y_pvt;
L1_postCal.ac_vel_z_pvt = rx_vel_z_pvt;

L1_postCal.ac_vel_x = rx_vel_xyz(:,1);
L1_postCal.ac_vel_y = rx_vel_xyz(:,2);
L1_postCal.ac_vel_z = rx_vel_xyz(:,3);

L1_postCal.ac_roll_pvt = rx_roll_deg_pvt;
L1_postCal.ac_pitch_pvt = rx_pitch_deg_pvt;
L1_postCal.ac_yaw_pvt = rx_yaw_deg_pvt;

L1_postCal.ac_roll = rx_attitude(:,1);
L1_postCal.ac_pitch = rx_attitude(:,2);
L1_postCal.ac_yaw = rx_attitude(:,3);

L1_postCal.rx_clk_bias_pvt = rx_clk_bias_m_pvt;
L1_postCal.rx_clk_drift_pvt = rx_clk_drift_mps_pvt;

L1_postCal.rx_clk_bias = rx_clk(:,1);
L1_postCal.rx_clk_drift = rx_clk(:,2);

% Part 1 ends

%% Part 2: Derive TX related variables
% 1 - derive TX positions and velocities
% 2 - map between GPS PRN and SVN
clc

% initalise variables 
tx_pos_x = zeros(J,I)+invalid;      tx_vel_x = zeros(J,I)+invalid;
tx_pos_y = zeros(J,I)+invalid;      tx_vel_y = zeros(J,I)+invalid;
tx_pos_z = zeros(J,I)+invalid;      tx_vel_z = zeros(J,I)+invalid;

tx_clk_bias = zeros(I)+invalid;

prn_code = zeros(J,I)+invalid;      sv_num = zeros(J,I)+invalid;

for i = 1:I

    % gps timestamp
    ddm_gps_timestamp1.gps_week = gps_week(i); 
    ddm_gps_timestamp1.gps_tow = gps_tow(i);

    for j = 1:J
    
        % PRN and SVN
        prn1 = transmitter_id(j,i);
        sv_num1 = SV_PRN_LUT(SV_PRN_LUT(:,1)==prn1,2);

        if ~isnan(prn1) && (prn1 ~= 0)
            [tx_pos_xyz1,tx_vel_xyz1,tx_clk_bias1,~] = gps_posvel(prn1,ddm_gps_timestamp1, ...
                gps_orbit_filename);

            tx_pos_x(j,i) = tx_pos_xyz1(1); tx_vel_x(j,i) = tx_vel_xyz1(1);
            tx_pos_y(j,i) = tx_pos_xyz1(2); tx_vel_y(j,i) = tx_vel_xyz1(2);
            tx_pos_z(j,i) = tx_pos_xyz1(3); tx_vel_z(j,i) = tx_vel_xyz1(3);

            tx_clk_bias(j,i) = tx_clk_bias1;

            prn_code(j,i) = prn1;           sv_num(j,i) = sv_num1;

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

% Part 2 ends

%% Part 3: L1a calibration
% 1- convert from raw counts to signal power in watts
clc

% TODO: 
% Calibration factors are determined by assuming the NGRx is constant
% and temperature independent
% Confusion on the small raw-count values

% initialise variables for L1a results
power_analog = zeros(5,40,J,I)+invalid;
ddm_ant = zeros(J,I)+invalid;
%inst_noise = zeros(I,1)+invalid;

for i = 1:I

    % retrieve noise standard deviation in counts for all three channels
    noise_std1 = [noise_std_dev_rf1(i),noise_std_dev_rf2(i),noise_std_dev_rf3(i)];

    for j = 1:J

        rf_source1 = rf_source(j,i);

        if ~isnan(rf_source1)

            ANZ_port1 = get_ANZ_port(rf_source1);
            raw_counts1 = raw_counts(:,:,j,i);                  % raw nadir antenna measurement in counts           

            % perform L1a calibration
            [signal_power_watts1,noise_power_watts1,snr_db1] = L1a_counts2watts(raw_counts1,ANZ_port1,noise_std1);
            inst_gain1 = max(raw_counts1,[],'all')/max(signal_power_watts1,[],'all');

            power_analog(:,:,j,i) = signal_power_watts1;
            ddm_ant(j,i) = ANZ_port1;

        end
        
    end    
end

% save outputs to L1 structure
L1_postCal.power_analog = power_analog;
L1_postCal.ddm_ant = ddm_ant;
%L1_postCal(i).noise_floor = mean(mean(noise_power_watts1));
%L1_postCal(i).ddm_snr = snr_db1;
%L1_postCal(i).inst_gain = inst_gain1;

% Part 3 ends

%% Part 4: L1b processing
clc

% initialise variables
sx_pos_x = zeros(J,I)+invalid;
sx_pos_y = zeros(J,I)+invalid;
sx_pos_z = zeros(J,I)+invalid;

sx_lat = zeros(J,I)+invalid;
sx_lon = zeros(J,I)+invalid;
sx_alt = zeros(J,I)+invalid;

sx_inc_angle = zeros(J,I)+invalid;
sx_dis_to_coast_km = zeros(J,I)+invalid;

tx_to_sp_range = zeros(J,I)+invalid;
rx_to_sp_range = zeros(I,1)+invalid;

gps_boresight = zeros(J,I)+invalid;

sx_theta_body = zeros(J,I)+invalid;
sx_az_body = zeros(J,I)+invalid;

sx_theta_enu = zeros(J,I)+invalid;
sx_az_enu = zeros(J,I)+invalid;

gps_tx_power_db_w = zeros(J,I)+invalid;
gps_ant_gain_db_i = zeros(J,I)+invalid;
static_gps_eirp = zeros(J,I)+invalid;

sx_rx_gain = zeros(J,I)+invalid;
cross_pol = zeros(J,I)+invalid;

brcs_ddm_sp_bin_delay_row = zeros(J,I)+invalid;
brcs_ddm_sp_bin_dopp_col = zeros(J,I)+invalid;

brcs = zeros(5,40,J,I)+invalid;
eff_scatter = zeros(5,40,J,I)+invalid;
nbrcs_scatter_area = zeros(J,I)+invalid;
nbrcs = zeros(J,I)+invalid;

les_scatter_area = zeros(J,I)+invalid;
ddm_les = zeros(J,I)+invalid;

tes_scatter_area = zeros(J,I)+invalid;
ddm_tes = zeros(J,I)+invalid;

surface_reflectivity = zeros(5,40,J,I)+invalid;
surface_reflectivity_peak = zeros(J,I)+invalid;

coherent_ratio = zeros(J,I)+invalid;
coherent_status = zeros(J,I)+invalid;

fresnel_minor = zeros(5,40,J,I)+invalid;
fresnel_major = zeros(5,40,J,I)+invalid;
fresnel_orientation = zeros(5,40,J,I)+invalid;

for i = 1:I

    % Part 4.0: retrieve variables for L1b processing
    % retrieve per-sample and write to structures
    rx1.rx_pos_xyz = rx_pos_xyz(i,:);
    rx1.rx_vel_xyz = rx_vel_xyz(i,:);
    
    rx1.rx_attitude = rx_attitude(i,:);
    rx1.rx_clk_drift = rx_clk(i,2);

    ddm1.delay_bin_res = delay_bin_res(i);
    ddm1.doppler_bin_res = doppler_bin_res(i);
    
    ddm1.delay_center_bin = delay_center_bin(i);
    ddm1.doppler_center_bin = doppler_center_bin(i);

    ddm1.num_delay_bins = num_delay_bins(i);
    ddm1.num_doppler_bins = num_doppler_bins(i);
    
    for j = 1:J

        % retrieve per_DDM variables
        % retrieve tx variables and write to tx structure
        tx_pos_xyz1 = [tx_pos_x(j,i) tx_pos_y(j,i) tx_pos_z(j,i)];
        tx_vel_xyz1 = [tx_vel_x(j,i) tx_vel_y(j,i) tx_vel_z(j,i)];
        prn_code1 = prn_code(j,i);
        sv_num1 = sv_num(j,i);        
        
        tx1.tx_pos_xyz = tx_pos_xyz1;
        tx1.tx_vel_xyz = tx_vel_xyz1;
        tx1.prn_code = prn_code1;
        tx1.sv_num = sv_num1;

        % retrieve ddm variables and write to ddm structure
        ddm1.raw_counts = raw_counts(:,:,j,i);

        ddm1.delay_center_chips = delay_center_chips(j,i);
        ddm1.doppler_center_Hz = doppler_center_Hz(j,i);

        ddm_ant1 = ddm_ant(j,i);

        % tracked direct delay chips, correct to within [0, 1023]
        add_range_to_sp_chips1 = meter2chips(add_range_to_sp(j,i));
        delay_dir_chips2 = delay_center_chips(j,i)-add_range_to_sp_chips1;
        delay_dir_chips1 = delay_correction(delay_dir_chips2,1023);

        ddm1.delay_dir_chips = delay_dir_chips1;
        ddm1.add_range_to_sp = add_range_to_sp_chips1;

        % analog power from L1a
        power_analog1 = power_analog(:,:,j,i);

        % only process these with valid prn code
        if prn_code1 ~= invalid

            % Part 4.1: SP solver
            % derive SP positions, angle of incidence and distance to coast
            % in kilometer
            [sx_pos_xyz1,inc_angle_deg1,dis_to_coast_km1] = sp_solver(tx1,rx1,ddm1, ...
                dem1,dem2,dtu10,landmask_nz,wgs84_flag);
            sx_pos_lla1 = ecef2lla(sx_pos_xyz1);            % <lat,lon,alt> of the specular reflection

            % define nadir antenna pattern to be used
            % installed antenna orientation angle is not included
            if ddm_ant1 == 2
                nadir_pattern = LHCP_pattern;

            elseif ddm_ant1 == 3
                nadir_pattern = RHCP_pattern;

            end

            % only process these with valid sx positions, i.e., in LOS
            if sum(sx_pos_xyz1) == invalid*3
                continue

            else
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
                    cross_pol1 = rx_rad1(1)-rx_rad1(2);

                elseif ddm_ant1 == 3
                    sx_rx_gain1 = rx_rad1(2);
                    cross_pol1 = rx_rad1(2)-rx_rad1(1);

                end

                % save to variables
                sx_pos_x(j,i) = sx_pos_xyz1(1);
                sx_pos_y(j,i) = sx_pos_xyz1(2);
                sx_pos_z(j,i) = sx_pos_xyz1(3);

                sx_inc_angle(j,i) = inc_angle_deg1;
                sx_dis_to_coast_km(j,i) = dis_to_coast_km1;

                sx_theta_body(j,i) = sx_angle_body1(1);
                sx_az_body(j,i) = sx_angle_body1(2);

                sx_theta_enu(j,i) = sx_angle_enu1(1);
                sx_az_enu(j,i) = sx_angle_enu1(2);

                gps_boresight = theta_gps1;

                tx_to_sp_range(j,i) = ranges1(1);
                rx_to_sp_range(j,i) = ranges1(2);

                gps_tx_power_db_w(j,i) = gps_rad1(1);
                gps_ant_gain_db_i(j,i) = gps_rad1(2);
                static_gps_eirp(j,i) = gps_rad1(3);

                sx_rx_gain(j,i) = sx_rx_gain1;                  % nadir antenna gain - main polarisation
                cross_pol(j,i) = cross_pol1;

                % Part 4.1 ends

                % Part 4.2: BRCS, NBRCS, and effective scattering area (A_eff)
                % derive BRCS
                brcs1 = ddm_brcs(power_analog1,gps_eirp_watt1, ...
                    sx_rx_gain1,R_tsx1,R_rsx1);

                % retrieve a discretised local surface for deriving A_eff
                L_local = 6090; res_local = 30;
                local_dem1 = get_local_dem(sx_pos_lla1,L_local,res_local, ...
                    dem1,dem2,dtu10,landmask_nz);

                T_coh = coherent_duration(i);                   % coherent duration

                % derive A_eff,nbrcs,LES and TES
                [A_eff1,sp_delay_bin_float1,sp_doppler_bin_float1] = ddm_Aeff(tx1,rx1,sx_pos_xyz1, ...
                    ddm1,local_dem1,T_coh);
                [nbrcs1,LES1,TES1] = ddm_nbrcs(brcs1,A_eff1,sp_delay_bin_float1,sp_doppler_bin_float1);

                % save to variables
                brcs_ddm_sp_bin_delay_row(j,i) = sp_delay_bin_float1-1;     % correct to 0-indexed
                brcs_ddm_sp_bin_dopp_col(j,i) = sp_doppler_bin_float1-1;    % correct to 0-indexed

                brcs(:,:,j,i) = brcs1;
                eff_scatter(:,:,j,i) = A_eff1;  
                nbrcs_scatter_area(j,i) = nbrcs1(1);
                nbrcs = nbrcs1(2);

                les_scatter_area(j,i) = LES1(1);
                ddm_les(j,i) = LES1(2);
                
                tes_scatter_area(j,i) = TES1(1);
                ddm_tes(j,i) = TES1(2);

                % Part 4.2 ends
                
                % Part 4.3: reflectivity and peak reflectivity
                [surf_refl1,surf_refl_peak1] = ddm_refl(power_analog1,gps_eirp_watt1, ...
                    sx_rx_gain1,R_tsx1,R_rsx1);

                % save to variables
                surface_reflectivity(:,:,j,i) = surf_refl1;
                surface_reflectivity_peak(j,i) = surf_refl_peak1;

                % Part 4.3 ends

                % Part 4.4: coherency detector
                [CR1,CS1] = coh_det(ddm1.raw_counts,snr_db);

                coherent_ratio(j,i) = CR1;
                coherent_status(j,i) = CS1;

                % Part 4.4 ends

                % Part 4.5: Fresnel-zone solver
                [fresnel_major1,fresnel_minor1,fresnel_orientation1] = fresnel(tx1,rx1,sx_pos_xyz1,inc_angle_deg1);

                % save to variables
                fresnel_major(j,i) = fresnel_major1;
                fresnel_minor(j,i) = fresnel_minor1;
                fresnel_orientation(j,i) = fresnel_orientation1;

                % Part 4.5 ends

            end

        end

    end    
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is part is for the variables that may need the value from a previous
% or a future timestamp, including the sx velocity and also quality flags
clc

for i = 1:I

    timestamp1 = L1_postCal(i).ddm_timestamp_utc;
    timestamp2 = L1_postCal(i+1).ddm_timestamp_utc;

    sx_pos_xyz1 = [ L1_postCal(i).sp_pos_x
                    L1_postCal(i).sp_pos_y
                    L1_postCal(i).sp_pos_z];

    sx_pos_xyz2 = [ L1_postCal(i+1).sp_pos_x
                    L1_postCal(i+1).sp_pos_y
                    L1_postCal(i+1).sp_pos_z];

    if (all(sx_pos_xyz1(:) == [-999,-999,-999])) || ...
            (all(sx_pos_xyz2(:) == [-999,-999,-999]))
        sx_vel_xyz = [-999 -999 -999];
    else
        sx_vel_xyz = (sx_pos_xyz2-sx_pos_xyz1)/(timestamp2-timestamp1);
    end

    L1_postCal(i).sx_vel_x = sx_vel_xyz(1);
    L1_postCal(i).sx_vel_y = sx_vel_xyz(2);
    L1_postCal(i).sx_vel_z = sx_vel_xyz(3);

    % intialise quality flag1
    quality_flag1 = zeros(29,1);

    % small and large sc attitude error
    roll = abs(L1_postCal(i).ac_roll);
    pitch = abs(L1_postCal(i).ac_pitch);
    yaw = abs(L1_postCal(i).ac_yaw);

    if (roll >= 30) || (pitch >= 10) || (yaw >= 5)
        quality_flag1(4) = 1;
    else
        quality_flag1(3) = 1;
    end

    % low_confidence_ddm_noise_floor and large_step_noise_floor
    noise_floor1 = L1_postCal(i).noise_floor;
    noise_floor2 = L1_postCal(i-1).noise_floor;

    noise_diff = abs(noise_floor1-noise_floor2)/noise_floor1;
    noise_diff_db = abs(pow2db(noise_floor1)-pow2db(noise_floor2));

    if noise_diff > 0.1
        quality_flag1(9) = 1;
    end

    if noise_diff_db > 0.24
        quality_flag1(12) = 1;
    end
    
    % sp over land and near land
    dis_to_coast = L1_postCal(i).dis_to_coast;

    if dis_to_coast >= 0
        quality_flag1(10) = 1;
    end

    if dis_to_coast >= -25
        quality_flag1(11) = 1;
    end

    % direct signal in ddm
    delay_dir_chips = L1_postCal(i).add_range_to_sp;
    delay_center_chips1 = delay_center_chips(i);

    d = delay_dir_chips-delay_center_chips;

    if abs(d) <= 4
        quality_flag1(14) = 1;
    end

    % sp_bin_delay_error
    sp_delay_bin_float = L1_postCal(i).brcs_ddm_sp_bin_delay_row;
    if (sp_delay_bin_float < 15) || (sp_delay_bin_float > 25)
        quality_flag1(17) = 1;
    end    
    
    % sp_bin_doppler_error
    sp_doppler_bin_float = L1_postCal(i).brcs_ddm_sp_bin_doppler_col;
    if (sp_doppler_bin_float < 2) || (sp_doppler_bin_float > 4)
        quality_flag1(18) = 1;
    end    
    
    % neg_ddma_brcs
    brcs = L1_postCal(i).brcs;
    ddma_brcs = brcs(floor(sp_doppler_bin_float)-1:floor(sp_doppler_bin_float)+1, ...
        floor(sp_delay_bin_float):floor(sp_delay_bin_float)+4);
    index = find(ddma_brcs<0);

    if ~isempty(index)
        quality_flag1(19) = 1;
    end

    % gps_pvt_sp3_error
    tx_pos_xyz = [  L1_postCal(i).tx_pos_x
                    L1_postCal(i).tx_pos_y
                    L1_postCal(i).tx_pos_z];

    tx_vel_xyz = [  L1_postCal(i).tx_vel_x
                    L1_postCal(i).tx_vel_y
                    L1_postCal(i).tx_vel_z];

    if (~all(tx_pos_xyz,'all')) || (~all(tx_vel_xyz,'all'))
        quality_flag1(20) = 1;
    end
    
    % sp_non_existent_error
    if all(sx_pos_xyz1(:) == [-999,-999,-999])
        quality_flag1(21) = 1;
    end

    % always false for this field
    quality_flag1(26) = 1;
    
    % ac_altitdue_out_of_nominal_range
    rx_alt = L1_postCal(i).ac_alt;
    if (rx_alt < 0) || (rx_alt > 15e3)
        quality_flag1(27) = 1;
    end

    % overall quality flag
    if (quality_flag1(2) == 1 || ...
            quality_flag1(4) == 1 || ...
            quality_flag1(5) == 1 || ...
            quality_flag1(6) == 1 || ...
            quality_flag1(7) == 1 || ...
            quality_flag1(8) == 1 || ...
            quality_flag1(12) == 1 || ...
            quality_flag1(13) == 1 || ...
            quality_flag1(14) == 1 || ...
            quality_flag1(15) == 1 || ...
            quality_flag1(16) == 1 || ...
            quality_flag1(17) == 1 || ...
            quality_flag1(18) == 1 || ...
            quality_flag1(19) == 1 || ...
            quality_flag1(21) == 1 || ...
            quality_flag1(22) == 1 || ...
            quality_flag1(23) == 1 || ...
            quality_flag1(24) == 1 || ...
            quality_flag1(25) == 1 || ...
            quality_flag1(26) == 1 || ...
            quality_flag1(27) == 1 || ...
            quality_flag1(28) == 1 || ...
            quality_flag1(29) == 1)
        quality_flag1(1) = 1;

    end

end

% L1 calibration ends

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

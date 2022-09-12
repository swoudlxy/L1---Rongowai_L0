% This script operates the L1b calibration for Rongowai
% Version: Matlab - 1.0
% inital version verified by CYGNSS L1
% Version: Matlab - 1.1
% 1) updated after testing with Rongowai simulated scenario data
% 2) convert all external source files date files

%% Load L0 data
clear
clc

% load L0 netCDF
path1 = '../dat/raw/';
L0_filename = [path1 'sample_flight_3_with_gps_time_corrected_L0.nc'];

% retrieve L0 data
time_utc_s = double(ncread(L0_filename,'time_utc_s'));                                      % utc time epoch
pvt_gps_week = double(ncread(L0_filename,'/science/gps_week'));
pvt_gps_sec = double(ncread(L0_filename,'/science/gps_seconds'));

%gps_vel_timestamp = double(ncread(L0_filename,'/eng/gps_velocity_time_stamp'));
%gps_vel_timestamp = gps_vel_timestamp(gps_vel_timestamp>0);

% PVT GPS week and sec - confirmed no att time
% GPS timestamps have invalid values in L0 testing file
%pvt_gps_week = double(ncread(L0_filename,'/science/GPS_week_of_SC_attitude'));
%pvt_gps_sec = double(ncread(L0_filename,'/science/GPS_second_of_SC_attitude'));

% rx positions in ECEF, metres
rx_pos_x_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_position_x_ecef_m'));
rx_pos_y_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_position_y_ecef_m'));
rx_pos_z_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_position_z_ecef_m'));

rx_pos_xyz_pvt = [rx_pos_x_pvt,rx_pos_y_pvt,rx_pos_z_pvt];

% rx velocity in ECEF, m/s
% rx velocity have invalid values or 0 in L0 testing file
rx_vel_x_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_velocity_x_ecef_mps'));
rx_vel_y_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_velocity_y_ecef_mps'));
rx_vel_z_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_velocity_z_ecef_mps'));

rx_vel_xyz_pvt = [rx_vel_x_pvt,rx_vel_y_pvt,rx_vel_z_pvt];

% rx attitude, deg
rx_pitch_deg_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_attitude_pitch_deg'));
rx_roll_deg_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_attitude_roll_deg'));
rx_yaw_deg_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_attitude_yaw_deg'));

rx_attitude_pvt = [rx_roll_deg_pvt,rx_pitch_deg_pvt,rx_yaw_deg_pvt];

% rx clock bias and drifts
rx_clk_bias_m_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_clock_bias_m'));
rx_clk_drift_mps_pvt = double(ncread(L0_filename,'/geometry/receiver/rx_clock_drift_mps'));

rx_clk_pvt = [rx_clk_bias_m_pvt,rx_clk_drift_mps_pvt];

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

% NGRx estimate additional delay path
add_range_to_sp = double(ncread(L0_filename,'/science/ddm/additional_range_to_SP'));

%% define external data paths and filenames
clc

% define IGS orbits filename (*.sp3)
% note this path is defined in C++ format
gps_orbit_filename = '..\\dat\\orbits\\igr21526.sp3';

% load DTU10, SRTM land DEM, and ocean-land mask (distance to coast)
dtu_path = '../dat/dtu/';
dtu_filename = 'dtu10_v1.dat';

dtu10 = get_dtu10([dtu_path dtu_filename]);

dem_path = '../dat/dem/';
load([dem_path,'NZSRTM30_v1.mat'])

% landmask/dist_to_coast has been updated to .dat format
landmask_path = '../dat/cst/';
landmask_filename = 'dist_to_coast_nz_v1.dat';
dist_to_coast_nz = dist_to_coast_mask([landmask_path landmask_filename]);

% load PRN-SV and SV-EIRP(static) LUT
% both LUTs have been converted to .dat format
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

% define azimuth and elevation angle in the antenna frame
res = 0.1;                                          % resolution in degrees
azim_deg = 0:res:360; elev_deg = 120:-1*res:0;

%% L1a calibration/processing
clc

% TODO: 
% Calibration factors are determined by assuming the NGRx is constant
% and temperature independent.
% Converting to incidence power on the antenna is not performed as more 
% parameters are required (e.g., cable loss, active antenna pattern, 
% LNA-temperature dependent gain variation, thermal model, etc).
% Noise standard deviation still needs to be redefined according to the
% latest L0 data.

% initialise a strucuture to save L1a results
L1_postCal = struct;

for i = 1:I

    rf_source1 = rf_source(i);

    % retrieve DDM raw counts at each timestamp
    raw_counts2 = raw_counts(:,:,:,i);
    raw_counts1 = squeeze(raw_counts2);         % remove dimension = 1

    % retrieve noise standard deviation
    noise_std1 = [noise_bin_avg_rf1,noise_std_dev_rf2,noise_std_dev_rf3];

    % perform L1a calibration
    [signal_power_watts1,noise_power_watts1,snr_db1] = counts2watts(raw_counts1,rf_source1,noise_std1);
    inst_gain1 = max(raw_counts1,'all')/max(signal_power_watts1,'all');

    % save outputs to L1 structure
    L1_postCal(i).power_analog = signal_power_watts1;
    L1_postCal(i).noise_floor = mean(mean(noise_power_watts1));
    L1_postCal(i).ddm_snr = snr_db1;
    L1_postCal(i).inst_gain = inst_gain1;

end

% L1a calibration ends

%% From PVT to DDM time epoch
clc

I = length(pvt_gps_sec);                % total length of samples

% initialise data array for timestamps
pvt_utc = zeros(I,1);   ddm_utc = zeros(I,1);
gps_week = zeros(I,1);  gps_tow = zeros(I,1);

% initialise a structure to save L1b results
L1_postCal = struct;

% derive and save ddm_timestamp_utc/gps
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
    [h,m,s] = hms(D_ddm1);
    [ddm_gps_week1,ddm_gps_sec1] = utc2gpstime(year,month,day,h,m,s);

    % save pvt and ddm timestamps for interpolation
    pvt_utc(i) = pvt_utc1;          ddm_utc(i) = ddm_utc1;
    gps_week(i) = ddm_gps_week1;    gps_tow(i) = ddm_gps_sec1;

    % save to L1 structure
    L1_postCal(i).pvt_timestamp_gps_week = pvt_gps_week1;
    L1_postCal(i).pvt_timestamp_gps_sec = pvt_gps_sec1;
    L1_postCal(i).pvt_timestamp_utc = pvt_utc1;

    L1_postCal(i).ac_pos_x_pvt = rx_pos_x_pvt(i);
    L1_postCal(i).ac_pos_y_pvt = rx_pos_y_pvt(i);
    L1_postCal(i).ac_pos_z_pvt = rx_pos_z_pvt(i);

    L1_postCal(i).ac_vel_x_pvt = rx_vel_x_pvt(i);
    L1_postCal(i).ac_vel_x_pvt = rx_vel_y_pvt(i);
    L1_postCal(i).ac_vel_x_pvt = rx_vel_z_pvt(i);

    L1_postCal(i).ac_roll_pvt = rx_roll_deg_pvt(i);
    L1_postCal(i).ac_pitch_pvt = rx_pitch_deg_pvt(i);
    L1_postCal(i).ac_yaw_pvt = rx_yaw_deg_pvt(i);

    L1_postCal(i).rx_clk_bias_pvt = rx_clk_bias_m_pvt(i);
    L1_postCal(i).rx_clk_drift_pvt = rx_clk_drift_mps_pvt(i);

    L1_postCal(i).ddm_timestamp_gps_week = ddm_gps_week1;
    L1_postCal(i).ddm_timestamp_gps_sec = ddm_gps_sec1;
    L1_postCal(i).ddm_timestamp_utc = ddm_utc1;

    L1_postCal(i).ddm_pvt_bias = non_coherent1/2;

end

% linear interpolation all the values at ddm timestamp
[rx_pos_xyz,rx_vel_xyz,rx_attitude,rx_clk] = PVT2ddm_timestamp(pvt_utc,ddm_utc,rx_pos_xyz_pvt, ...
    rx_vel_xyz_pvt,rx_attitude_pvt,rx_clk_pvt);

rx_clk_bias_m = interp1(pvt_utc,rx_clk_bias_m_pvt,ddm_utc);
rx_clk_drift_mps = interp1(pvt_utc,rx_clk_drift_mps_pvt,ddm_utc);

% save global duration values
L1_postCal.time_coverage_start = datetime(ddm_utc(1),'ConvertFrom','posixtime');        % human time
L1_postCal.time_coverage_end = datetime(ddm_utc(end),'ConvertFrom','posixtime');
L1_postCal.time_coverage_duration = ddm_utc(end)-ddm_utc(1);

L1_postCal.dem_source = 'SRTM30';

%% L1b calibration/processing
clc

for i = 1:I

    % Part 0: retrieve all time-dependent parameters and variables for calibration
    
    % retrieve parameters at each timestamp
    % product of L1a
    power_analog1 = L1_postCal(i).power_analog;
    ddm_snr1 = L1_postCal(i).ddm_snr;

    % gps timestamp
    gps_timestamp1.gps_week = gps_week(i);
    gps_timestamp1.gps_tow = gps_tow(i);

    % transmitter ID, PRN
    tx_id1 = transmitter_id(i);
    
    % rx-related parameters to rx structure
    rx1.rx_pos_xyz = rx_pos_xyz(i,:);
    rx1.rx_vel_xyz = rx_vel_xyz(i,:);
    rx1.rx_attitude = rx_attitude(i,:);
    rx1.rx_clock_drift_mps = rx_clk_drift_mps(i);

    % rx <lat,lon,alt> coordinate
    rx_pos_lla1 = ecef2lla(rx_pos_xyz(i,:));
    
    % ddm-related parameters to ddm structure
    raw_counts2 = raw_counts(:,:,:,i);
    ddm1.raw_counts = squeeze(raw_counts2);         % remove dimensions = 1
 
    ddm1.delay_dir_chips = delay_dir_chips(i);
    
    ddm1.delay_bin_res = delay_bin_res(i);
    ddm1.doppler_bin_res = doppler_bin_res(i);

    ddm1.delay_center_bin = delay_center_bin(i);
    ddm1.doppler_center_bin = doppler_center_bin(i);

    ddm1.delay_center_chips = delay_center_chips(i);
    ddm1.doppler_ceneter_Hz = doppler_center_Hz(i);

    ddm1.num_delay_bins = num_delay_bins(i);
    ddm1.num_doppler_bins = num_doppler_bins(i);

    rf_source1 = rf_source(i);                      % rf3 = LHCP, rf5 = RHCP

    % derive GNSS satellite position and velocity at each timestamp
    [tx_pos_xyz1,tx_vel_xyz1,gps_clk_bias1,gps_clk_drift1] = gps_posvel(tx_id1, ...
        gps_timestamp1,gps_orbit_filename);

    % tx-related parameters to tx structure
    tx1.tx_pos_xyz = tx_pos_xyz1;
    tx1.tx_vel_xyz = tx_vel_xyz1;
    tx1.tx_id = tx_id1;

    % save rx drift and bias to L1b structure
    L1_postCal(i).rx_clk_bias = rx_clk_bias_m(i,1);
    L1_postCal(i).rx_clk_drift = rx_clk_drift_mps(i,2);

    % TODO:
    % rx_clk_bias_rate_pvt does not exist in the current L0
    % rx_clk_bias_rate needs to be interpolated based on that

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Part 1: SP solver
    % Part 1 computes the coordinate of specular reflections and other 
    % SP-related variables such as angles in variable frames, ranges,
    % and etc.
    clc

    % compute SP locations, incidence angle and distance to coast
    [sx_pos_xyz1,inc_angle_deg1,dis_to_coast_km1] = sp_solver(tx1,rx1,ddm1, ...
        dtu10,NZSRTM30,dist_to_coast_nz);

    % compute SP-related variables for each solved SP
    % define LHCP or RHCP pattern to use:
    % rf_source 4 = NGRx_RF3_NB = ANZ_RF2 = nadir_LHCP
    % rf_source 8 = NGRx_RF5_NB = ANZ_RF3 = nadir_RHCP
    if rf_source1 == 4
        nadir_pattern = LHCP_pattern;
    elseif rf_source1 == 8
        nadir_pattern = RHCP_pattern;
    end

    % initialise variables
    num_sx = size(sx_pos_xyz1,1);                   %number of SP positions solved by SP solver

    sp_theta_body = zeros(num_sx,1);    sp_az_body = zeros(num_sx,1);
    sp_theta_enu = zeros(num_sx,1);     sp_az_enu = zeros(num_sx,1);
    gps_off_boresight = zeros(num_sx,1);

    range_tsx = zeros(num_sx,1);        range_rsx = zeros(num_sx,1);

    gps_pow = zeros(num_sx,1);          gps_gain_dbi = zeros(num_sx,1);
    stat_eirp_watt = zeros(num_sx,1);

    rx_gain = zeros(num_sx,1);
    cross_pol = zeros(num_sx,1);

    for l = 1:num_sx

        sx_pos_xyz2 = sx_pos_xyz1(l,:);

        [sp_angle_body1,sp_angle_enu1,theta_gps1,range1,gps_rad1,rx_rad1] = spRelated(tx1,rx1,sx_pos_xyz2, ...
            SV_PRN_LUT,SV_eirp_LUT,nadir_pattern);

        sp_theta_body(l) = sp_angle_body1(1);   sp_az_body(l) = sp_angle_body1(2);
        sp_theta_enu(l) = sp_angle_enu1(1);     sp_az_enu(l) = sp_angle_enu1(2);
        gps_off_boresight(l) = theta_gps1;

        range_tsx(l) = range1(1);               range_rsx(l) = range1(2);

        gps_pow(l) = gps_rad1(1);               gps_gain_dbi = gps_rad1(2);
        stat_eirp_watt(l) = gps_rad1(3);

        lhcp_gain_dbi = rx_rad1(1);             rhcp_gain_dbi = rx_rad1(2);

        % compute cross polarisation
        if rf_source1 == 4
            rx_gain(l) = lhcp_gain_dbi;
            cross_pol(l) = lhcp_gain_dbi-rhcp_gain_dbi;

        elseif rf_source1 == 8
            rx_gain(l) = rhcp_gain_dbi;
            cross_pol(l) = rhcp_gain_dbi-lhcp_gain_dbi;

        end        

    end

    % save SP sovler outputs to the post-calibrated structure
    L1_postCal(i).tx_pos_x = tx_pos_xyz1(1);
    L1_postCal(i).tx_pos_y = tx_pos_xyz1(2);
    L1_postCal(i).tx_pos_z = tx_pos_xyz1(3);

    L1_postCal(i).tx_vel_x = tx_vel_xyz1(1);
    L1_postCal(i).tx_vel_y = tx_vel_xyz1(2);
    L1_postCal(i).tx_vel_z = tx_vel_xyz1(3);

    L1_postCal(i).ac_pos_x = rx_pos_xyz(i,1);
    L1_postCal(i).ac_pos_y = rx_pos_xyz(i,2);
    L1_postCal(i).ac_pos_z = rx_pos_xyz(i,3);

    L1_postCal(i).ac_lat = rx_pos_lla1(1);
    L1_postCal(i).ac_lon = rx_pos_lla1(2);
    L1_postCal(i).ac_alt = rx_pos_lla1(3);    

    L1_postCal(i).ac_roll = rx_attitude(i,1);
    L1_postCal(i).ac_pitch = rx_attitude(i,2);
    L1_postCal(i).ac_yaw = rx_attitude(i,3); 

    L1_postCal(i).sp_pos_x = sx_pos_xyz1(:,1);
    L1_postCal(i).sp_pos_y = sx_pos_xyz1(:,2);
    L1_postCal(i).sp_pos_z = sx_pos_xyz1(:,3);

    L1_postCal(i).sp_inc_angle = inc_angle_deg1;
    L1_postCal(i).dis_to_coast = dis_to_coast_km1;

    L1_postCal(i).sp_theta_body = sp_theta_body;
    L1_postCal(i).sp_az_body = sp_az_body;

    L1_postCal(i).sp_theta_enu = sp_theta_enu;
    L1_postCal(i).sp_az_enu = sp_az_enu;

    L1_postCal(i).tx_to_sp_range = range_tsx;
    L1_postCal(i).rx_to_sp_range = range_rsx;

    L1_postCal(i).gps_off_boresight_angle_deg = gps_off_boresight;
    L1_postCal(i).static_gps_eirp = stat_eirp_watt;
    L1_postCal(i).gps_tx_power_db_w = gps_pow;

    L1_postCal(i).ddm_source = rf_source1;
    L1_postCal(i).sp_rx_gain = rx_gain;
    L1_postCal(i).cross_pol = cross_pol;

    L1_postCal(i).zenith_code_phase = delay_dir_chips(i);
        
    % Part 1: SP solver ends

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Part 2: BRCS, A_eff, and NBRCS
    % Part 2 computes the bistatic radar cross section (BRCS), effective
    % scattering area (A_eff) and normalised BRCS
    clc

    % take average of the parameters used in the radar equation for
    % computation, these values are not reported in the L1b product
    gps_eirp_mean = mean(stat_eirp_watt);
    rx_gain_mean = mean(rx_gain);

    range_tsx_mean = mean(range_tsx);
    range_rsx_mean = mean(range_rsx);

    T_coh = coherent_duration(i);
    
    % compute BRCS
    brcs1 = ddm_brcs(power_analog1,gps_girp_mean,rx_gain_mean,range_tsx_mean,range_rsx_mean);

    % compute A_eff
    num_grid = 201;                     % number of grids for each side for the local dem

    % define local dem according to the solved SP coordinate(s)
    if num_sx == 1
        sx_pos_xyz_local = sx_pos_xyz1;

    elseif num_sx > 1
        sx_pos_xyz_local(1) = mean(sx_pos_xyz1(:,1));
        sx_pos_xyz_local(2) = mean(sx_pos_xyz1(:,2));
        sx_pos_xyz_local(3) = mean(sx_pos_xyz1(:,3));

    end

    local_dem = localdem(sx_pos_xyz_local,num_grid,NZSRTM30,dtu10,dist_to_coast_nz);
    [A_eff1,sp_delay_bin_float1,sp_doppler_bin_float1] = ddm_Aeff(tx1,rx1,sx_pos_xyz1, ...
        ddm1,local_dem,T_coh);

    % compute DDMA NBRCS
    [nbrcs1,LES1,TES1] = ddm_nbrcs(brcs1,A_eff1,sp_delay_bin_float1,sp_doppler_bin_float1);

    % save outputs to the post-calibration structure
    L1_postCal(i).brcs = brcs1;
    L1_postCal(i).eff_scatter = A_eff1;

    L1_postCal(i).brcs_ddm_sp_bin_delay_row = sp_delay_bin_float1;
    L1_postCal(i).brcs_ddm_sp_bin_dopp_col = sp_doppler_bin_float1;

    L1_postCal(i).ddm_nbrcs = nbrcs1.nbrcs_value;
    L1_postCal(i).nbrcs_scatter_area = nbrcs1.nbrcs_scatter;

    L1_postCal(i).LES_scatter = LES1.LES_scatter;
    L1_postCal(i).LES_slope = LES1.LES_slope;

    L1_postCal(i).TES_scatter = TES1.TES_scatter;
    L1_postCal(i).TES_slope = TES1.TES_slope;

    % Part 2: BRCS, A_eff, and NBRCS ends

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Part 3: Reflectivity and peak reflectivity

    % Part 3 derives the surface reflectivity and the peak surface
    % reflectivity. This computation is straightforward and using Friis
    % transmission equation
    clc

    [refl1,peak_refl1] = ddm_refl(power_analog1,gps_girp_mean,rx_gain_mean, ...
        range_tsx_mean,range_rsx_mean);

    % save outputs to post-calibration structure
    L1_postCal(i).reflectivity = refl1;
    L1_postCal(i).peak_reflectivity = peak_refl1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Part 4: Coherent detection and Fresnel dimension

    % Part 4 derives coherent ratio and state based on the input raw counts
    % Fresnel dimensions are also derived here for coherency ratio >= 2
    clc

    % derive coherency ratio (CR) and coherency state(CS)
    raw_counts1 = squeeze(raw_counts2);
    [CR1,CS1] = coh_det(raw_counts1,snr_db1);

    % derive Fresnel dimensions for CR >= 2
    % it is reasonable to assume only one SP exist for coherent reflections
    if CR1 >= 2
        [major_axis1,minor_axis1,orientation1] = fresnel(tx1,rx1,sx_pos_xyz1,inc_angle_deg1);

    elseif CR1 < 2
        major_axis1 = -999;
        minor_axis1 = -999;
        orientation1 = -999;
        
    end

    % save output to post-calibration structure
    L1_postCal(i).coherency_ratio = CR1;
    L1_postCal(i).coherency_state = CS1;
    
    L1_postCal(i).fresnel_major = major_axis1;
    L1_postCal(i).fresnel_minor = minor_axis1;
    L1_postCal(i).fresnel_direction = orientation1;

    % Part 4: Coherent detection and Fresnel dimension ends

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Part 5: Diagnostic values
    clc
    
    add_range_to_sp_chips = meter2chips(add_range_to_sp(i));
    L1_postCal(i).add_range_to_sp = add_range_to_sp_chips;

    % Part 5 ends

end

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

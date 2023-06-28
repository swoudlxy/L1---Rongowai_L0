% version B of L1 code to process mulitple L0 files at a time
clear
clc

% load L0 data
L0_path = '../dat/raw/';
L0_filenames = dir([L0_path '*.nc']);

L = length(L0_filenames);

%% load external files
clc

% load L1a calibration tables
L1a_path = '../dat/L1a_cal/';

L1a_cal_ddm_counts_db = readmatrix([L1a_path 'L1A_cal_ddm_counts_dB_v1.dat']);  % version control - 28 June
L1a_cal_ddm_power_dbm = readmatrix([L1a_path 'L1A_cal_ddm_power_dBm_v1.dat']);

% load SRTM DEM
dem_path = '../dat/dem/';
dem_file = 'dem_200m.dat';

dem = get_dem2([dem_path dem_file]); 

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
% update antenna pattern LUT - 28 June
rng_path = '../dat/rng/';

LHCP_L_filename = 'GNSS_LHCP_L_gain_db_i_v2.dat';
LHCP_R_filename = 'GNSS_LHCP_R_gain_db_i_v2.dat';

RHCP_L_filename = 'GNSS_RHCP_L_gain_db_i_v2.dat';
RHCP_R_filename = 'GNSS_RHCP_R_gain_db_i_v2.dat';

[~,~,LHCP_L_gain_db_i] = get_ant_pattern([rng_path,LHCP_L_filename]);
[~,~,LHCP_R_gain_db_i] = get_ant_pattern([rng_path,LHCP_R_filename]);
LHCP_pattern.LHCP = LHCP_L_gain_db_i;
LHCP_pattern.RHCP = LHCP_R_gain_db_i;

[~,~,RHCP_L_gain_db_i] = get_ant_pattern([rng_path,RHCP_L_filename]);
[~,~,RHCP_R_gain_db_i] = get_ant_pattern([rng_path,RHCP_R_filename]);
RHCP_pattern.LHCP = RHCP_L_gain_db_i;
RHCP_pattern.RHCP = RHCP_R_gain_db_i;

% scattering area LUT
% LUT change - 28 June
A_phy_LUT_files = dir('../dat/A_phy_LUT/*.dat');
path = A_phy_LUT_files(1).folder;
L = length(A_phy_LUT_files);

% read input variables
[rx_alt_bins,inc_angle_bins,az_angle_bins] = get_A_phy_inputs([path '/' A_phy_LUT_files(12).name]);

num_rx_alt = length(rx_alt_bins);
num_inc_angle = length(inc_angle_bins);
num_az_angle = length(az_angle_bins);

A_phy_LUT_all = struct;

for l = 1:L-1

    A_phy_LUT_filename = A_phy_LUT_files(l).name;
    A_phy_LUT_file = [path '/' A_phy_LUT_filename];
    A_phy_LUT = get_A_phy_LUT(A_phy_LUT_file,num_rx_alt,num_inc_angle,num_az_angle);

    sp_doppler_frac1 = str2double(A_phy_LUT_filename(11:end-4));
    sp_doppler_frac_range = -0.5:0.1:0.5;

    % ordering from smallest to largest fraction
    [~,k] = min(abs(sp_doppler_frac_range-sp_doppler_frac1));

    A_phy_LUT_all(k).sp_doppler_frac = sp_doppler_frac1;
    A_phy_LUT_all(k).A_phy_LUT = A_phy_LUT;

end

% L1 dictionary name
L1_dict_name = '../dat/L1_Dict/L1_Dict_v2_2.xlsx';

%% get post-calibrated L1 product
clc

for l = 2

    filename = L0_filenames(l).name;
    path = L0_filenames(l).folder;

    L0_filename = [path '/' filename];

    L1_postCal = get_L1_product(L0_filename, ...
                                L1a_cal_ddm_counts_db,L1a_cal_ddm_power_dbm, ...
                                dem,dtu10,landmask_nz,lcv_mask,water_mask, ...
                                SV_PRN_LUT,SV_eirp_LUT, ...
                                LHCP_pattern,RHCP_pattern, ...
                                rx_alt_bins,inc_angle_bins,az_angle_bins,A_phy_LUT_all);

    % the below saves solved L1 products as a MATLAB structure, not necessary for daily processing  
    %save(['../out/L1_postCalData/' filename(1:end-3) '_L1.mat'],'L1_postCal','-v7.3');

    % the below packets the L1 products as a netCDF
    L1_netCDF = ['../out/L1_netCDFs/' filename(1:end-3) '_L1.nc'];
    sample_info = get_netcdf(L1_netCDF,L1_dict_name,L1_postCal);
    
    disp(l);    % display the number of files being processed

end
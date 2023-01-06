% version B of L1 code to process mulitple L0 files
clear
clc

% load L0 data
L0_path = '../dat/raw/Week1_Nov/';
L0_filenames = dir([L0_path '*.nc']);

L = length(L0_filenames);

%% load external files
clc

% load L1a calibration tables
L1a_path = '../dat/L1a_cal/';

L1a_cal_ddm_counts_db = readmatrix([L1a_path 'L1A_cal_ddm_counts_dB.dat']);
L1a_cal_ddm_power_dbm = readmatrix([L1a_path 'L1A_cal_ddm_power_dBm.dat']);

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

% physical element size LUT
phy_ele_size = readmatrix('../dat/dem/phy_ele_size.dat');

%% get post-calibrated L1 product
clc

for l = 24:L

    filename = L0_filenames(l).name;
    path = L0_filenames(l).folder;

    L0_filename = [path '/' filename];

    L1_postCal = get_L1_product(L0_filename, ...
                                L1a_cal_ddm_counts_db,L1a_cal_ddm_power_dbm, ...
                                dem,dtu10,landmask_nz,lcv_mask,water_mask, ...
                                SV_PRN_LUT,SV_eirp_LUT, ...
                                LHCP_pattern,RHCP_pattern, ...
                                phy_ele_size);

    save(['../out/L1_postCalData/' filename(1:end-3) '_L1.mat'],'L1_postCal','-v7.3');

end











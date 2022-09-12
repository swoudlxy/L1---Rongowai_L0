%clear workspace
close all
clear
clc

%define paths
path1 = '../dat/raw/';
path2 = '../int/';

file_name = 'cyg01.ddmi.s20201231-000000-e20201231-234659.l1.land-reflectivity.sand032.nc';
file_name_ocean = 'cyg01.ddmi.s20201231-000000-e20201231-235959.l1.power-brcs.a30.d31.nc';

filename = [path1 file_name];
filename_ocean = [path1 file_name_ocean];

%below are retrieved from land sandbox
ddm_gps_sec = double(ncread(filename,'ddm_timestamp_gps_sec'));

tx_pos_x = double(ncread(filename,'tx_pos_x')); 
tx_pos_y = double(ncread(filename,'tx_pos_y')); 
tx_pos_z = double(ncread(filename,'tx_pos_z')); 

tx_vel_x = double(ncread(filename,'tx_vel_x'));
tx_vel_y = double(ncread(filename,'tx_vel_y'));
tx_vel_z = double(ncread(filename,'tx_vel_z'));

sc_pos_x = double(ncread(filename,'sc_pos_x')); 
sc_pos_y = double(ncread(filename,'sc_pos_y')); 
sc_pos_z = double(ncread(filename,'sc_pos_z'));

sc_vel_x = double(ncread(filename,'sc_vel_x'));
sc_vel_y = double(ncread(filename,'sc_vel_y'));
sc_vel_z = double(ncread(filename,'sc_vel_z'));

%sp parameters for reference
sp_pos_x = double(ncread(filename,'sp_pos_x')); 
sp_pos_y = double(ncread(filename,'sp_pos_y'));
sp_pos_z = double(ncread(filename,'sp_pos_z'));

sp_lat = double(ncread(filename,'sp_lat')); 
sp_lon = double(ncread(filename,'sp_lon'));
sp_alt = double(ncread(filename,'sp_alt'));

sp_vel_x = double(ncread(filename,'sp_vel_x'));
sp_vel_y = double(ncread(filename,'sp_vel_y'));
sp_vel_z = double(ncread(filename,'sp_vel_z'));

%additional code phase delay
add_range_to_sp_L1 = double(ncread(filename,'add_range_to_sp'));

%zenith code phase
zenith_code_phase_L1 = double(ncread(filename,'zenith_code_phase'));

%absolute code phase and delay at compressed DDM center
delay_center_L1 = double(ncread(filename,'code_phase_comp_ddm_center_row'));
dopp_center_L1 = double(ncread(filename,'doppler_comp_ddm_center_col'));

%doppler shifts due to the instrument
rx_dopp_L1 = double(ncread(filename,'rx_clk_doppler'));

%flight software delay and doppler shifts
fsw_comp_delay_shift_L1 = double(ncread(filename,'fsw_comp_delay_shift'));
fsw_comp_dopp_shift_L1 = double(ncread(filename,'fsw_comp_dopp_shift'));

%raw counts
raw_counts_L1 = double(ncread(filename,'raw_counts'));

%below are retrieved from public ocean L1 product
ddm_gps_sec_ocean = double(ncread(filename_ocean,'ddm_timestamp_gps_sec'));

sp_pos_x_ocean = double(ncread(filename_ocean,'sp_pos_x')); 
sp_pos_y_ocean = double(ncread(filename_ocean,'sp_pos_y'));
sp_pos_z_ocean = double(ncread(filename_ocean,'sp_pos_z'));

sp_inc_angle_ocean = double(ncread(filename_ocean,'sp_pos_z'));

%% data for Channel 2
clc

%tx rx coordinates
tx_xyz_ch2 = [tx_pos_x(2,:)'    tx_pos_y(2,:)'  tx_pos_z(2,:)'];
rx_xyz_ch2 = [sc_pos_x          sc_pos_y        sc_pos_z];

%sp coordinates for reference
sp_xyz_ch2 = [sp_pos_x(2,:)'    sp_pos_y(2,:)'  sp_pos_z(2,:)'];

%tx rx sp velocity
tx_vel_xyz_ch2 = [tx_vel_x(2,:)'    tx_vel_y(2,:)'  tx_vel_z(2,:)'];
rx_vel_xyz_ch2 = [sc_vel_x          sc_vel_y        sc_vel_z];
sp_vel_xyz_ch2 = [sp_vel_x(2,:)'    sp_vel_y(2,:)'  sp_vel_z(2,:)'];

add_range_to_sp_ch2 = add_range_to_sp_L1(2,:)';

zenith_code_phase_ch2 = zenith_code_phase_L1(2,:)';

delay_center_ch2 = delay_center_L1(2,:)';
dopp_center_ch2 = dopp_center_L1(2,:)';

rx_dopp_ch2 = rx_dopp_L1(2,:)';

raw_counts_ch2 = raw_counts_L1(:,:,2,:);
raw_counts_ch2 = squeeze(raw_counts_ch2);

%remove invalid rows 
index = find(~isnan(tx_xyz_ch2(:,1)));
%index = find(sp_lat(2,:)>=13 & sp_lat(2,:)<=15 & sp_lon(2,:)>=32 & sp_lon(2,:)<=34);

tx_xyz = tx_xyz_ch2(index,:);
rx_xyz = rx_xyz_ch2(index,:);
sp_xyz = sp_xyz_ch2(index,:);

tx_vel_xyz = tx_vel_xyz_ch2(index,:);
rx_vel_xyz = rx_vel_xyz_ch2(index,:);

add_range_to_sp = add_range_to_sp_ch2(index);

zenith_code_phase = zenith_code_phase_ch2(index);

delay_center = delay_center_ch2(index);
dopp_center = dopp_center_ch2(index);

rx_dopp = rx_dopp_ch2(index);

raw_counts = raw_counts_ch2(:,:,index);

%retrieve same timestamp data from public ocean netCDF
ddm_gps_sec_land = ddm_gps_sec(index);

I = length(ddm_gps_sec_land);
sp_xyz_oceanL1 = zeros(I,3);

for i = 1:I
    
    ddm_gps_sec_land1 = ddm_gps_sec_land(i);
    index2 = find(ddm_gps_sec_ocean == ddm_gps_sec_land1);

    sp_pos_x1 = sp_pos_x_ocean(2,index2);
    sp_pos_y1 = sp_pos_y_ocean(2,index2);
    sp_pos_z1 = sp_pos_z_ocean(2,index2);

    sp_xyz_oceanL1(i,:) = [sp_pos_x1 sp_pos_y1 sp_pos_z1];

end



%%
save([path2 'tx_xyz.mat'],'tx_xyz');
save([path2 'rx_xyz.mat'],'rx_xyz');
save([path2 'sp_xyz.mat'],'sp_xyz');

save([path2 'tx_vel_xyz.mat'],'tx_vel_xyz');
save([path2 'rx_vel_xyz.mat'],'rx_vel_xyz');

save([path2 'zenith_code_phase.mat'],'zenith_code_phase');

save([path2 'delay_center.mat'],'delay_center');
save([path2 'dopp_center.mat'],'dopp_center');

save([path2 'rx_dopp.mat'],'rx_dopp');

save([path2 'raw_counts.mat'],'raw_counts');
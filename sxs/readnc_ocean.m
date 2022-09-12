%clear workspace
close all
clear
clc

%define paths
path1 = '../dat/raw/';
path2 = '../int/';

file_name = 'cyg01.ddmi.s20201231-000000-e20201231-235959.l1.power-brcs.a30.d31.nc';
filename = [path1 file_name];

%gps weeks and secs
gps_week = double(ncread(filename,'ddm_timestamp_gps_week'));
gps_sec = double(ncread(filename,'ddm_timestamp_gps_sec'));

%tx and sc positions and velocity
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

%sc altitude
sc_roll = double(ncread(filename,'sc_roll'));
sc_pitch = double(ncread(filename,'sc_pitch'));
sc_yaw = double(ncread(filename,'sc_yaw'));

%sp positions
sp_pos_x = double(ncread(filename,'sp_pos_x')); 
sp_pos_y = double(ncread(filename,'sp_pos_y'));
sp_pos_z = double(ncread(filename,'sp_pos_z'));

prn_code_L1 = double(ncread(filename,'prn_code'));

sp_inc_angle_L1 = double(ncread(filename,'sp_inc_angle'));         %incidence angle

gps_off_boresight_angle_deg_L1 = double(ncread(filename,'gps_off_boresight_angle_deg'));

sp_theta_orbit_L1 = double(ncread(filename,'sp_theta_orbit'));
sp_az_orbit_L1 = double(ncread(filename,'sp_az_orbit'));

sp_theta_body_L1 = double(ncread(filename,'sp_theta_body'));
sp_az_body_L1 = double(ncread(filename,'sp_theta_body'));

static_gps_eirp_L1 = double(ncread(filename,'static_gps_eirp'));

%% data for Channel 2
clc

%tx rx coordinates
tx_pos_xyz = [tx_pos_x(2,:)'    tx_pos_y(2,:)'  tx_pos_z(2,:)'];
rx_pos_xyz = [sc_pos_x          sc_pos_y        sc_pos_z];

sc_att     = [sc_roll           sc_pitch        sc_yaw];

%tx rx sp velocity
tx_vel_xyz = [tx_vel_x(2,:)'    tx_vel_y(2,:)'  tx_vel_z(2,:)'];
rx_vel_xyz = [sc_vel_x          sc_vel_y        sc_vel_z];

%sp coordinates for reference
sp_pos_xyz = [sp_pos_x(2,:)'    sp_pos_y(2,:)'  sp_pos_z(2,:)'];

prn_code = prn_code_L1(2,:)';

sp_inc_angle = sp_inc_angle_L1(2,:)';

gps_off_boresight_angle_deg = gps_off_boresight_angle_deg_L1(2,:)';

sp_theta_orbit = sp_theta_orbit_L1(2,:)';
sp_az_orbit = sp_az_orbit_L1(2,:)';

sp_theta_body = sp_theta_body_L1(2,:)';
sp_az_body = sp_az_body_L1(2,:)';

static_gps_eirp = static_gps_eirp_L1(2,:)';

%% save matrix
clc

save([path2 'tx_pos_xyz.mat'],'tx_pos_xyz');
save([path2 'rx_pos_xyz.mat'],'rx_pos_xyz');

save([path2 'sc_att.mat'],'sc_att');

save([path2 'sp_pos_xyz.mat'],'sp_pos_xyz');

save([path2 'tx_vel_xyz.mat'],'tx_vel_xyz');
save([path2 'rx_vel_xyz.mat'],'rx_vel_xyz');

save([path2 'sp_inc_angle.mat'],'sp_inc_angle');

save([path2 'gps_off_boresight_angle_deg.mat'],'gps_off_boresight_angle_deg')

save([path2 'sp_theta_orbit.mat'],'sp_theta_orbit');
save([path2 'sp_az_orbit.mat'],'sp_az_orbit');

save([path2 'sp_theta_body.mat'],'sp_theta_body');
save([path2 'sp_az_body.mat'],'sp_az_body');

save([path2 'static_gps_eirp.mat'],'static_gps_eirp')

%% save as txt to feed C++
clc

txt_path = '../GPS_SP3_STG/';
fop = fopen([txt_path,'inputs.txt'],'w');

input_matrix = [prn_code,ddm_timestamp_gps_week,ddm_timestamp_gps_sec];
[M,N] = size(input_matrix);

for m = 1:M
    for n = 1:N

        data = input_matrix(m,n);
        fprintf(fop,'%s ',num2str(data));
        
    end
    fprintf(fop,'\n');
end

fclose(fop);
























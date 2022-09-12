%clear workspace
close all
clear
clc

%define paths
path1 = '../dat/raw/';
path2 = '../int/';

file_name = '00001_1546304357_NZWN_NZGS_RLK488-1546307905-c819de.nc4';
filename = [path1 file_name];

%tx and sc positions and velocity
tx_pos_x = double(ncread(filename,'/geometry/transmitter/tx_position_x_ecef_m')); 
tx_pos_y = double(ncread(filename,'/geometry/transmitter/tx_position_y_ecef_m')); 
tx_pos_z = double(ncread(filename,'/geometry/transmitter/tx_position_z_ecef_m')); 

sc_pos_x = double(ncread(filename,'/geometry/receiver/rx_position_x_ecef_m')); 
sc_pos_y = double(ncread(filename,'/geometry/receiver/rx_position_y_ecef_m')); 
sc_pos_z = double(ncread(filename,'/geometry/receiver/rx_position_z_ecef_m'));

%% data for Channel 4
clc

ch = 4;

%tx rx coordinates
tx_pos_xyz = [tx_pos_x(ch,:)'    tx_pos_y(ch,:)'  tx_pos_z(ch,:)'];
rx_pos_xyz = [sc_pos_x           sc_pos_y         sc_pos_z];


%% save matrix
clc

save([path2 'tx_pos_xyz.mat'],'tx_pos_xyz');
save([path2 'rx_pos_xyz.mat'],'rx_pos_xyz');
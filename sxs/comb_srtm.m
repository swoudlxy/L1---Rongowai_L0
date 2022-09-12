% this function combines four separate DEMs to one DEM

%srtm = comb_srtm(file_path,file1,file2,file3,file4)
clear
clc

file_path = '../dat/dem/';

% read the 1st file
file1 = 'nzsrtm_30_part1_v1.dat';
file2 = 'nzsrtm_30_part2_v1.dat';
file3 = 'nzsrtm_30_part3_v1.dat';
file4 = 'nzsrtm_30_part4_v1.dat';

srtm_dem1 = process_srtm(file_path,file1);
srtm_dem2 = process_srtm(file_path,file2);
srtm_dem3 = process_srtm(file_path,file3);
srtm_dem4 = process_srtm(file_path,file4);

lat1 = srtm_dem1.lat;   lat2 = srtm_dem3.lat;
lon1 = srtm_dem1.lon;   lon2 = srtm_dem2.lon; 

ele11 = srtm_dem1.ele;  ele12 = srtm_dem2.ele;
ele21 = srtm_dem3.ele;  ele22 = srtm_dem4.ele;

%%

srtm_lat = [lat1(1:end-1) lat2];
srtm_lon = [lon1(1:end-1) lon2];

srtm_ele = [ele11(1:end-1,1:end-1)  ele12(1:end-1,:);
            ele21(:,1:end-1)        ele22(:,:)];
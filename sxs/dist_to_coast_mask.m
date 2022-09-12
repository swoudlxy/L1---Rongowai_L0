% this function read distance to coast data files to the workspace
% input: path and filename of the dist_to_coast data file
% output: dist_to_coast_mask structure

function dist_to_coast = dist_to_coast_mask(filename)

fid = fopen(filename,'r');

lat_min = fread(fid,1,'double');
lat_max = fread(fid,1,'double');
num_lat = fread(fid,1,'uint16');

lon_min = fread(fid,1,'double');
lon_max = fread(fid,1,'double');
num_lon = fread(fid,1,'uint16');

num_pixel = num_lat*num_lon;
map_data = fread(fid,num_pixel,'double');

fclose(fid);

dist_to_coast.lat = linspace(lat_min,lat_max,num_lat);
dist_to_coast.lon = linspace(lon_min,lon_max,num_lon);
dist_to_coast.dist = reshape(map_data,num_lat,num_lon);
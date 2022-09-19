% this function outputs the DEM data

function srtm_dem = get_dem(filename)

fid = fopen(filename,'r');

lat_min = fread(fid,1,'double');
lat_max = fread(fid,1,'double');
num_lat = fread(fid,1,'uint16');

lon_min = fread(fid,1,'double');
lon_max = fread(fid,1,'double');
num_lon = fread(fid,1,'uint16');

num_pixel = num_lat*num_lon;
map_data = fread(fid,num_pixel,'uint16');

fclose(fid);

lat = linspace(lat_min,lat_max,num_lat);
lon = linspace(lon_min,lon_max,num_lon);
srtm_data = reshape(map_data,num_lat,num_lon);

srtm_dem.lat = lat;
srtm_dem.lon = lon;
srtm_dem.ele = srtm_data;
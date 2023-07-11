function dem = get_dem2(filename)

fid = fopen(filename);

lat_min = fread(fid,1,'float32');
lat_res = fread(fid,1,'float32');
num_lat = fread(fid,1,'float32');

lon_min = fread(fid,1,'float32');
lon_res = fread(fid,1,'float32');
num_lon = fread(fid,1,'float32');

num_grid = num_lat*num_lon;

ele = fread(fid,num_grid,'uint16');

fclose(fid);

lat = lat_min+((1:1:num_lat)-1)*lat_res;
lon = lon_min+((1:1:num_lon)-1)*lon_res;

ele = reshape(ele,num_lat,[]);

dem.lat = lat;
dem.lon = lon;
dem.ele = ele;
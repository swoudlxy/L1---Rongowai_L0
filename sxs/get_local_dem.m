% this function outputs the local DEM data around the local coordinate P
% P - LLA coordinate

function local_dem = get_local_dem(P,L,dem_res,dem_data1,dem_data2)

lat_P = P(1);   lon_P = P(2);

num_pixel = L/dem_res;
half_num_pixel = floor(num_pixel/2);

% sparse dem structures
lat1 = dem_data1.lat;   lon1 = dem_data1.lon;   ele1 = dem_data1.ele;
lat2 = dem_data2.lat;   lon2 = dem_data2.lon;   ele2 = dem_data2.ele;

% construct a local elevation map from srtm_30
if (lat_P>=lat1(end-half_num_pixel)) 
    [~,lat_index] = min(abs(lat1-lat_P));
    [~,lon_index] = min(abs(lon1-lon_P));

    local_lat = lat1(lat_index-half_num_pixel:lat_index+half_num_pixel);
    local_lon = lon1(lon_index-half_num_pixel:lon_index+half_num_pixel);
    local_ele = ele1(lat_index-half_num_pixel:lat_index+half_num_pixel, ...
        lon_index-half_num_pixel:lon_index+half_num_pixel);
    
elseif (lat_P<=lat2(half_num_pixel))
    [~,lat_index] = min(abs(lat2-lat_P));
    [~,lon_index] = min(abs(lon2-lon_P));

    local_lat = lat2(lat_index-half_num_pixel:lat_index+half_num_pixel);
    local_lon = lon2(lon_index-half_num_pixel:lon_index+half_num_pixel);
    local_ele = ele2(lat_index-half_num_pixel:lat_index+half_num_pixel, ...
        lon_index-half_num_pixel:lon_index+half_num_pixel);

elseif (lat_P<lat1(end-half_num_pixel)) && (lat_P>=lat1(end))
    [~,lat_index] = min(abs(lat1-lat_P));
    [~,lon_index] = min(abs(lon1-lon_P));

    diff = num_pixel-length(lat1(lat_index-half_num_pixel:end));
    
    local_lat = [lat1(lat_index-half_num_pixel:end) lat2(2:diff+1)];
    local_lon = lon1(lon_index-half_num_pixel:lon_index+half_num_pixel);
    local_ele = [ele1(lat_index-half_num_pixel:end,lon_index-half_num_pixel:lon_index+half_num_pixel);
        ele2(2:diff+1,lon_index-half_num_pixel:lon_index+half_num_pixel)];

elseif (lat_P>lat2(half_num_pixel)) && (lat_P<=lat2(1))
    [~,lat_index] = min(abs(lat2-lat_P));
    [~,lon_index] = min(abs(lon2-lon_P));

    diff = num_pixel-length(lat2(1:lat_index+half_num_pixel));
    
    local_lat = [lat1(end-diff:end-1) lat2(1:lat_index+half_num_pixel)];
    local_lon = lon2(lon_index-half_num_pixel:lon_index+half_num_pixel);
    local_ele = [ele1(end-diff:end-1,lon_index-half_num_pixel:lon_index+half_num_pixel);
        ele2(1:lat_index+half_num_pixel,lon_index-half_num_pixel:lon_index+half_num_pixel)];

end

local_dem.lat = local_lat;
local_dem.lon = local_lon;
local_dem.ele = local_ele;

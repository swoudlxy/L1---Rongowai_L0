% this function outputs the local DEM data around the local coordinate P
% P - LLA coordinate

function local_dem = get_local_dem(P,L,res,dem_data,dtu_model,dist_to_coast)

ocean_land_margin = 0;
lat_P = P(1);   lon_P = P(2);

num_pixels = L/res;
half_num_pixel = floor(num_pixels/2);

% sparse dem structures
lat = dem_data.lat;
lon = dem_data.lon;
ele = dem_data.ele;

[~,lat_index] = min(abs(lat-lat_P));
[~,lon_index] = min(abs(lon-lon_P));

%step = lon(2)-lon(1);
%lat_index = ceil((lat(1)-lat_P)/step);
%lon_index = ceil((lon_P-lon(1))/step);

local_lat = lat(lat_index-half_num_pixel:lat_index+half_num_pixel);
local_lon = lon(lon_index-half_num_pixel:lon_index+half_num_pixel);

if dist_to_coast > ocean_land_margin
    local_ele = ele(lat_index-half_num_pixel:lat_index+half_num_pixel, ...
        lon_index-half_num_pixel:lon_index+half_num_pixel);

else
    local_ele = zeros(num_pixels);
    
    for i = 1:num_pixels
        for j = 1:num_pixels

            pixel_lat = local_lat(i);
            pixel_lon = local_lon(j);
            pixel_ele = get_map_value(pixel_lat,pixel_lon,dtu_model);

            local_ele(i,j) = pixel_ele;

        end
    end

end

local_dem.lat = local_lat;
local_dem.lon = local_lon;
local_dem.ele = double(local_ele);
% this function returns the mss of a pixel <lat1, lon1>
% The mss of the pixel is derived by interpolating a DTU10 datum
% with 1 deg resolution
% Inputs:
% 1) <lat,lon>: geo coordinate of the pixel to be computed
% 2) model: dtu model, currently using DTU10
% Output:
% 1) mss: linearly interpolated mean sea surface elevation

function ele = get_map_value(lat,lon,model)

lat_map = model.lat;    lat_res = abs(lat_map(2)-lat_map(1));
lon_map = model.lon;    lon_res = abs(lon_map(2)-lon_map(1));
ele_map = model.ele;

y0 = lat;  
x0 = lon;

% longitude adjustment
if lon < 0
    lon = lon+360;
elseif lon > 360
    lon = lon-360;
end

% get four corners coordinates and elevations
y1_index = ceil(abs((lat-lat_map(1)))/lat_res); y1 = lat_map(y1_index);
y2_index = y1_index+1;

x1_index = ceil(abs((lon-lon_map(1)))/lon_res); x1 = lon_map(x1_index);
x2_index = x1_index+1;

ele1 = ele_map(y1_index,x1_index);
ele2 = ele_map(y1_index,x2_index);
ele3 = ele_map(y2_index,x1_index);
ele4 = ele_map(y2_index,x2_index);

% interpolation
fy = y0-y1;
fx = x0-x1;

temp1 = ele1*(1-fx)+fx*ele2;
temp2 = ele3*(1-fx)+fx*ele4;
ele = temp1*(1-fy)+temp2*fy;
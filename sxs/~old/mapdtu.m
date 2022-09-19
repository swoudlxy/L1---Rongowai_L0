% this function returns the mss of a pixel <lat1, lon1>
% The mss of the pixel is derived by interpolating a DTU10 datum
% with 1 deg resolution
% Inputs:
% 1) <lat,lon>: geo coordinate of the pixel to be computed
% 2) model: dtu model, currently using DTU10
% Output:
% 1) mss: linearly interpolated mean sea surface elevation

function mss = mapdtu(lat,lon,model)

lat_dtu = model.lat;
lon_dtu = model.lon;
mss_dtu = model.mss;

%longitude adjustment
if lon < 0
    lon = lon+360;
elseif lon > 360
    lon = lon-360;
end

y0 = lat;
x0 = lon;

%identify four close pixels
x1 = floor(x0); x2 = ceil(x0);
y1 = floor(y0); y2 = ceil(y0);

x1_index = find(lon_dtu == x1); x2_index = find(lon_dtu == x2);
y1_index = find(lat_dtu == y1); y2_index = find(lat_dtu == y2);

mss1 = mss_dtu(y1_index,x1_index);
mss2 = mss_dtu(y1_index,x2_index);
mss3 = mss_dtu(y2_index,x1_index);
mss4 = mss_dtu(y2_index,x2_index);

fy = y0-y1;
fx = x0-x1;

temp1 = mss1*(1-fx)+fx*mss2;
temp2 = mss3*(1-fx)+fx*mss4;
mss = temp1*(1-fy)+temp2*fy;
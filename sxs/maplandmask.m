% this function returns the distance to coast of a pixel <lat1, lon1>
% The mss of the pixel is derived by interpolating a local NZ
% ocean-land mask with 1 minute resolution
% Inputs:
% 1) <lat,lon>: coordinate of the pixel to be computed
% 2) model: ocean-land mask model

function dist = maplandmask(lat,lon,model)

lat_landmask = model.lat;
lon_landmask = model.lon;
dist_landmask = model.dist;

y0 = lat;
x0 = lon;

%identify four close pixels
x1_index = find(lon_landmask<x0,1,'last'); x2_index = find(lon_landmask>x0,1,'first');
y1_index = find(lat_landmask<y0,1,'last'); y2_index = find(lat_landmask>y0,1,'first');

x1 = lon_landmask(x1_index);    x2 = lon_landmask(x2_index);    x = [x1,x2];
y1 = lat_landmask(y1_index);    y2 = lat_landmask(y2_index);    y = [y1;y2];

dist1 = dist_landmask(y1_index,x1_index);
dist2 = dist_landmask(y1_index,x2_index);
dist3 = dist_landmask(y2_index,x1_index);
dist4 = dist_landmask(y2_index,x2_index);

V = [dist1 dist2;dist3 dist4];

dist = interp2(x,y,V,x0,y0);      %interpolated result
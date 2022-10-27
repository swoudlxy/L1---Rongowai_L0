function dem = get_dem(filename)

[A,R] = readgeoraster(filename);

latlim = R.LatitudeLimits;
lonlim = R.LongitudeLimits;
rastersize = R.RasterSize;

lat = linspace(latlim(2),latlim(1),rastersize(1));
lon = linspace(lonlim(1),lonlim(2),rastersize(2));

dem.lat = lat;
dem.lon = lon;
dem.ele = A;
function pek_value = get_pek_value(lat,lon,water_mask)

pek = water_mask.pek;
ref = water_mask.ref;

grid_size = ref.RasterSize(1);
res_deg = 10/grid_size;

lat_max = ref.LatitudeLimits(2);
lon_min = ref.LongitudeLimits(1);

%lat_pek = linspace(ref.LatitudeLimits(1),ref.LatitudeLimits(2),grid_size);
%lon_pek = linspace(ref.LongitudeLimits(1),ref.LongitudeLimits(2),grid_size);

lat_index = ceil((lat_max-lat)/res_deg);
lon_index = ceil((lon-lon_min)/res_deg);

pek_value = pek(lat_index,lon_index);
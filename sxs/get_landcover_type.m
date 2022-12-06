% this function returns the landcover type of the coordinate P (lat lon)
% over landsurface

function landcover_type = get_landcover_type(lat,lon,lcv)

[M,N,~] = size(lcv);

lat_max = -34;      lat_range = 13.5;   lat_res = lat_range/M;
lon_min = 165.75;   lon_range = 13.5;   lon_res = lon_range/N;

lat_index = ceil((lat_max-lat)/lat_res);
lon_index = ceil((lon-lon_min)/lon_res);

lcv_RGB1 = double(lcv(lat_index,lon_index,:));
lcv_RGB = squeeze(lcv_RGB1)/255;

% RGB triplet
color = [   0.8     0       0.8;        % 1: artifical
            0.6     0.4     0.2;        % 2: barely vegetated
            0       0       1;          % 3: inland water
            1       1       0;          % 4: crop
            0       1       0;          % 5: grass
            0.6     0.2     0;          % 6: shrub
            0       0.2     0];         % 7: forest

% compare with the RGB triplet to infer the landcover type
I = size(color,1);
pixel_diff = zeros(I,1);

if sum(lcv_RGB) == 3
    landcover_type = -2;                % no data

elseif sum(lcv_RGB) < 3
        
    for i = 1:I
    
        color1 = color(i,:);
        diff1 = sum(abs(lcv_RGB'-color1));
        pixel_diff(i) = diff1;

    end

    temp = find(pixel_diff == min(pixel_diff));
    landcover_type = temp;

end
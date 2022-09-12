% this function reconstructs a local DEM centred at the coordinate of the
% input point. The reconstruct DEM reconstructs a DEM from both DTU 
% and SRTM database, i.e., the region identified as "land" uses SRTM 
% elevation values, whereas the region identified as "ocean" uses DTU10 
% mean sea surface (MSS) values 
% Inputs:
% 1) S_pos_xyz: ECEF position of the point S
% 2) num_grid: number of grids for one side
% 3) DTU, SRTM, landmask: DTU, SRTM and landmask structure
% Outputs:
% 1) local_dem: local reconstructed DEM structure, including local lat,
% local lon, and local elevation, centred at S

function local_dem = localdem(S_pos_xyz,num_grid,SRTM_model,DTU_model,landmask_model)

% sparse SRTM model structures
lat_srtm = SRTM_model.lat;
lon_srtm = SRTM_model.lon;
ele_srtm = SRTM_model.ele;

% reconstruct local dem
S_lla = ecef2lla(S_pos_xyz);
[~,S_lat_index] = min(abs(lat_srtm-S_lla(1)));
[~,S_lon_index] = min(abs(lon_srtm-S_lla(2)));

% <lat lon> range of local dem
% take one more pixel here to compute size
half_length = floor(num_grid/2)+1;

lat_local_range = S_lat_index-half_length:S_lat_index+half_length;
lon_local_range = S_lon_index-half_length:S_lon_index+half_length;

lat_local = lat_srtm(lat_local_range);
lon_local = lon_srtm(lon_local_range);
ele_local = ele_srtm(lat_local_range,lon_local_range);

% loop over local dem and replace ocean elevation with MSS
for m = 1:num_grid
    for n = 1:num_grid

        lat_grid = lat_local(m);
        lon_grid = lon_local(n);

        % replace the elevations identified as ocean with DTU10 MSS
        dist = maplandmask(lat_grid,lon_grid,landmask_model);

        if dist <= -25
            ele_grid = mapdtu(lat_grid,lon_grid,DTU_model);

        else
            continue

        end

        ele_local(m,n) = ele_grid;

    end
end

% loop over to derive physical size for each pixel
dA = zeros(num_grid);

for m = 2:num_grid-1
    for n = 2:num_grid-1
        
        point1 = [lat_local(m),lon_local(n)];
        
        % coordinates of four points around point1
        point2 = [lat_local(m-1),lon_local(n)];
        point3 = [lat_local(m),lon_local(n-1)];
        point4 = [lat_local(m+1),lon_local(n)];
        point5 = [lat_local(m),lon_local(n+1)];

        sidelength1 = m_lldist([point1(2) point2(2)],[point1(1) point2(1)]) ... 
            + m_lldist([point1(2) point4(2)],[point1(1) point4(1)]);
        sidelength2 = m_lldist([point1(2) point3(2)],[point1(1) point3(1)]) ... 
            + m_lldist([point1(2) point5(2)],[point1(1) point5(1)]);

        dA(m,n) = sidelength1*sidelength2/4;

    end
end



local_dem.lat = lat_local(2:end-1);
local_dem.lon = lon_local(2:end-1);
local_dem.ele = ele_local(2:end-1,2:end-1);
local_dem.dA = dA;




































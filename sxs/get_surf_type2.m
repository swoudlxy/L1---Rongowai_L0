% this function returns the surface type of a coordinate P <lat lon>

function surface_type = get_surf_type2(P,cst_mask,pek_mask,lcv_mask)

% sparse inland water mask
water_mask_160E_40S = pek_mask.water_mask_160E_40S;
water_mask_170E_30S = pek_mask.water_mask_170E_30S;
water_mask_170E_40S = pek_mask.water_mask_170E_40S;

P_lla = ecef2lla(P);
lat_P = P_lla(1); lon_P = P_lla(2);

landcover_type = get_landcover_type2(lat_P,lon_P,lcv_mask);

% determine if coordinate over inland water
lat_pek = floor(abs(lat_P/10))*10;
lon_pek = floor(abs(lon_P/10))*10;

file_id = [num2str(lon_pek) 'E_' num2str(lat_pek) 'S'];
water_mask1 = eval(['water_mask_' file_id]);
    
pek_value = get_pek_value(lat_P,lon_P,water_mask1);

dist_coast = get_map_value(lat_P,lon_P,cst_mask);

if (pek_value>0) && (landcover_type~=-1) && dist_coast>0.5
    surface_type = 3;            % coordinate on inland water
elseif (pek_value>0) && (dist_coast<0.5)
    surface_type = -1;
else
    surface_type = landcover_type;
end
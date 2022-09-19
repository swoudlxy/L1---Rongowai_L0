% SP solver derives the coordinate(s) of the specular reflection (sx)
% SP solver also reports the local incidence angle and the distance to coast in km where the SP occur
% All variables are reported in ECEF
% Inputs:
% 1) tx, rx, and ddm properties as a structure
% 2) DEM models: dtu10, NZSRTM30, and land-ocean mask
% 3) wgs84_flag = 1 controls the solved SP coordinates only on WGS84 datum
% Outputs:
% 1) sx_pos_xyz: sx position and velocity vectors in ECEF
% 2) in_angle_deg: local incidence angle at the specular reflection

function [sx_pos_xyz,inc_angle_deg,dis_to_coast_km] = sp_solver(tx,rx,ddm,dem1,dem2,dtu10,dist_to_coast_nz,wgs84_flag)

% retrieve tx and rx position vectors
tx_pos_xyz = tx.tx_pos_xyz;
rx_pos_xyz = rx.rx_pos_xyz;

% step 0 - check if LOS exists
LOS_flag = los_status(tx_pos_xyz,rx_pos_xyz);       % LOS flag, 1 = LOS exit

if LOS_flag == 1
        
    % step 1 - derive inital SP coordinate on WGS84
    [sx_xyz_coarse,sx_lla_coarse] = coarsetune(tx_pos_xyz,rx_pos_xyz);

    if wgs84_flag == 1
        sx_pos_xyz = sx_xyz_coarse;
        rsx = rx_pos_xyz-sx_pos_xyz;
        inc_angle_deg = acosd(dot(sx_pos_xyz,rsx)/(norm(sx_pos_xyz)*norm(rsx)));

    elseif wgs84_flag == 0
        % step 2 - land-ocean mask to determine landcover type
        dist = get_map_value(sx_lla_coarse(1),sx_lla_coarse(2),dist_to_coast_nz);

        % ocean-land boundary: -25 km
        ocean_land_margin = -25;
        if dist < ocean_land_margin
            land_flag = 0;                          % raw SP on ocean

        elseif dist >= ocean_land_margin
            land_flag = 1;                          % raw SP on land

        end

        % step 3 - fine-tune initial SP coordinate
        if land_flag == 0
    
            % step 3.1 - ocean SP
            L_ocean_deg = 1;                        % search area starting from 1*1 deg
            res_ocean_meter = 0.01;                 % targeted grid resolution - 1 cm

            [sx_pos_xyz,inc_angle_deg] = finetune_ocean(tx_pos_xyz,rx_pos_xyz, ...
                sx_lla_coarse,dtu10,L_ocean_deg,res_ocean_meter);
        
        elseif land_flag == 1

            L_land_meter = 6030;                    % search area boundary in meters
            res_land_meter = 30;                    % resolution of the DEM in meters

            [sx_pos_xyz,inc_angle_deg] = finetune_land(tx,rx,sx_lla_coarse, ...
                L_land_meter,res_land_meter,ddm,dem1,dem2);

        end

    end

    % step 4 - derive distance to coast in kilometer
    sx_lla = ecef2lla(sx_pos_xyz);    
    dis_to_coast_km = get_map_value(sx_lla(1),sx_lla(2),dist_to_coast_nz);

elseif LOS_flag == 0

    %no sx if no LOS between rx and tx
    sx_pos_xyz = [-9999,-9999,-9999];
    inc_angle_deg = -9999;
    dis_to_coast_km = -9999;

end
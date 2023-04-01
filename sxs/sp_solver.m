% SP solver derives the coordinate(s) of the specular reflection (sx)
% SP solver also reports the local incidence angle and the distance to coast in km where the SP occur
% All variables are reported in ECEF
% Inputs:
% 1) tx and rx positions
% 2) DEM models: dtu10, NZSRTM30, and land-ocean mask
% Outputs:
% 1) sx_pos_xyz: sx positions in ECEF
% 2) in_angle_deg: local incidence angle at the specular reflection
% 3) distance to coast in kilometer
% 4) LOS flag

function [sx_pos_xyz,inc_angle_deg,d_snell_deg,dist_to_coast_km,LOS_flag] = sp_solver(tx_pos_xyz,rx_pos_xyz,dem_data,dtu10,dist_to_coast_nz)

% step 0 - check if LOS exists
LOS_flag = los_status(tx_pos_xyz,rx_pos_xyz);

if LOS_flag == 1
        
    % step 1 - derive SP coordinate on WGS84 and DTU10
    [sx_xyz_coarse,sx_lla_coarse] = coarsetune(tx_pos_xyz,rx_pos_xyz);

    L_ocean_deg = 1;                            % initial searching region in degrees
    res_ocean_meter = 0.01;                     % converge criteria 0.01 meter

    sx_pos_xyz = finetune_ocean(tx_pos_xyz,rx_pos_xyz,sx_lla_coarse, ...
        dtu10,L_ocean_deg,res_ocean_meter);
    
    % derive local angles
    sx_pos_lla = ecef2lla(sx_pos_xyz);
    dist = get_map_value(sx_pos_lla(1),sx_pos_lla(2),dist_to_coast_nz);

    local_dem = get_local_dem(sx_pos_lla,90,30,dem_data,dtu10,dist);
    [theta_i,theta_s,phi_i,phi_s] = angles(local_dem.lat,local_dem.lon,local_dem.ele,tx_pos_xyz,rx_pos_xyz);
        
    % step 2 - project to local DEM for land SPs
    if dist > 0

        local_height = local_dem.ele;
        local_height = local_height(2,2);       % local height of the SP

        % projection to local dem
        term1 = sx_xyz_coarse/norm(sx_xyz_coarse);
        term2 = term1*local_height;
        sx_pos_xyz = sx_xyz_coarse+term2;

    end

    %inc_angle_deg = theta_i;
    
    v_tsx = tx_pos_xyz-sx_pos_xyz;
    unit_tsx = v_tsx/norm(v_tsx); 
    unit_sx = sx_pos_xyz/norm(sx_pos_xyz);
    inc_angle_deg = acosd(dot(unit_tsx,unit_sx));

    d_theta = theta_i-theta_s;
    d_phi1 = sind(phi_s-(phi_i+180))/cosd(phi_s-(phi_i+180));
    d_phi = atand(d_phi1);

    d_snell_deg = abs(d_theta)+abs(d_phi); 
    
    dist_to_coast_km = dist;

elseif LOS_flag == 0

    %no sx if no LOS between rx and tx
    sx_pos_xyz = [NaN,NaN,NaN];
    inc_angle_deg = NaN;
    d_snell_deg = NaN;
    dist_to_coast_km = NaN;

end
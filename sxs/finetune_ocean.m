% This function fine tunes the SP coordiantes using a DTU10 datum
% Inputs:
% 1) TX and 2) RX coordinates in the form of ECEF-XYZ and
% 3) SP coordinate in the form of LLA
% 4) model: earth model - currently DTU10
% 5) L: inital searching area in deg
% 6) res_grid: targeted resolution of each grid when quitting the iteration
% in metres
% Output: return 
% 1) fine-tuned SP coordiantes in ECEF-XYZ, and
% 2) local incidence angle

function [sx_xyz_final,theta_i] = finetune_ocean(tx_pos_xyz,rx_pos_xyz,sp_lla_coarse,model,L,res_grid)

%derive SP on the ocean surface
res = 1000;

while res > res_grid
    
    [res,~,sp_temp] = finetune(tx_pos_xyz,rx_pos_xyz,sp_lla_coarse,L,model);
    
    %parameters for the next iteration
    L = L*2/11;                                         %new searching area
    sp_lla_coarse = sp_temp;                            %new SP coordinate

end

sx_lla_final = sp_temp;                                 %finalised sp in lla
sx_xyz_final = lla2ecef(sx_lla_final);                  %finalised sp in ecef-xyz

%derive incidence angle
rsx = rx_pos_xyz-sx_xyz_final;
theta_i = acosd(dot(rsx,sx_xyz_final)/(norm(rsx)*norm(sx_xyz_final)));

%{
wgs84 = wgs84Ellipsoid;

%convert Tx coordinate to ENU coordinate centred at the computed sx
[tx_e,tx_n,tx_u] = ecef2enu(tx_pos_xyz(1),tx_pos_xyz(2),tx_pos_xyz(3), ...
    sx_lla_final(1),sx_lla_final(2),sx_lla_final(3),wgs84);
tx_enu = [tx_e,tx_n,tx_u];

unit_z = [0,0,1];                                       %unit vector along normal to the surface

%incidence angles at the SP
theta_i = acosd(dot(tx_enu,unit_z)/(norm(tx_enu)*norm(unit_z)));
%}
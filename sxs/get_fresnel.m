% this function derives Fresnel dimensions based on the Tx, Rx and Sx
% positions.
% Fresnel dimension is computed only the DDM is classified as coherent
% reflection.

function [fresnel_coeff,fresnel_axis,fresnel_orientation] = get_fresnel(tx_pos_xyz,...
    rx_pos_xyz,sx_pos_xyz,dist_to_coast,inc_angle,ddm_ant)

wgs84 = wgs84Ellipsoid('meter');
eps_ocean = 74.62+1j*51.92;         % complex permittivity of ocean

% sparse structures
%tx_pos_xyz = tx.tx_pos_xyz;
%rx_pos_xyz = rx.rx_pos_xyz;
%sx_pos_xyz = sx.sx_pos_xyz;
%dist_to_coast = sx.dist_to_coast;

% define constants
fc = 1575.42e6;             % operating frequency
c = 299792458;              % speed of light
lambda = c/fc;              % wavelength

% compute dimensions
R_tsp = norm(tx_pos_xyz-sx_pos_xyz);
R_rsp = norm(rx_pos_xyz-sx_pos_xyz);

term1 = R_tsp*R_rsp;
term2 = R_tsp+R_rsp;

% semi axis
a = sqrt(lambda*term1/term2);   % major semi
b = a/cosd(inc_angle);          % minor semi

% compute orientation relative to North
sx_lla = ecef2lla(sx_pos_xyz);

[tx_e,tx_n,~] = ecef2enu(tx_pos_xyz(1),tx_pos_xyz(2),tx_pos_xyz(3),sx_lla(1),sx_lla(2),sx_lla(3),wgs84);
[rx_e,rx_n,~] = ecef2enu(rx_pos_xyz(1),rx_pos_xyz(2),rx_pos_xyz(3),sx_lla(1),sx_lla(2),sx_lla(3),wgs84); 

tx_en = [tx_e,tx_n];
rx_en = [rx_e,rx_n];

vector_tr = rx_en-tx_en;
unit_north = [0,1];

term3 = dot(vector_tr,unit_north);
term4 = norm(vector_tr)*norm(unit_north);

theta = acosd(term3/term4);

fresnel_axis = [2*a 2*b];
fresnel_orientation = theta;

% fresenel coefficient only compute for ocean SPs
if dist_to_coast <= 0

    sint = sind(inc_angle);
    cost = cosd(inc_angle);

    temp1 = sqrt(eps_ocean-sint*sint);

    R_vv = (eps_ocean*cost-temp1)/(eps_ocean*cost+temp1);
    R_hh = (cost-temp1)/(cost+temp1);

    R_rl = (R_vv-R_hh)/2;
    R_rr = (R_vv+R_hh)/2;

    if ddm_ant == 2
        fresnel_coeff = abs(R_rl)*abs(R_rl);

    elseif ddm_ant == 3
        fresnel_coeff = abs(R_rr)*abs(R_rr);
    
    end

else
    fresnel_coeff = nan;

end














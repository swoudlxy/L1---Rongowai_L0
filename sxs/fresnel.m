% this function derives Fresnel dimensions based on the Tx, Rx and Sx
% positions.
% Fresnel dimension is computed only the DDM is classified as coherent
% reflection.

function [fresnel_major,fresnel_minor,fresnel_orientation] = fresnel(tx,rx,sx_pos_xyz,inc_angle)

wgs84 = wgs84Ellipsoid('meter');

% sparse structures
tx_pos_xyz = tx.tx_pos_xyz;
rx_pos_xyz = rx.rx_pos_xyz;

% define constants
fc = 1575.42e6;             % operating frequency
c = 299792458;              % speed of light
lambda = c/fc;              % wavelength

% compute dimensions
R_tsp = norm(tx_pos_xyz-sx_pos_xyz);
R_rsp = norm(rx_pos_xyz-sx_pos_xyz);

term1 = R_tsp*R_rsp;
term2 = R_tsp+R_rsp;

a = sqrt(lambda*term1/term2);
b = a/cosd(inc_angle);

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

fresnel_major = 2*a;
fresnel_minor = 2*b;
fresnel_orientation = theta;













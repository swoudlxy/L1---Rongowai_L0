% This function computes the local incidence and reflection angles of 
% the middle pixel in a 3 by 3 DEM pixel matrix
% Inputs: 
% 1) lat,lon, and ele matrices of the 3*3 pixel matrix
% 2) Tx and Rx coordinates ECEF(x,y,z)
% Outputs: 
% 1) theta_i, phi_i: local incidence angle along elevation and azimuth
% angles in degree
% 2) theta_s, phi_s: local scattering (reflection) angles along elevation
% and azimuth angles in degree

function [theta_i,theta_s,phi_i,phi_s] = angles(lat,lon,ele,tx_pos_xyz,rx_pos_xyz)

wgs84 = wgs84Ellipsoid;

s0 = [lat(2) lon(2) ele(2,2)];            %origin of the local enu frame
s0_enu = [0,0,0];

%convert tx and rx to local ENU centred at s0
[tx_e,tx_n,tx_u] = ecef2enu(tx_pos_xyz(1),tx_pos_xyz(2),tx_pos_xyz(3),s0(1),s0(2),s0(3),wgs84);
[rx_e,rx_n,rx_u] = ecef2enu(rx_pos_xyz(1),rx_pos_xyz(2),rx_pos_xyz(3),s0(1),s0(2),s0(3),wgs84);

tx_enu = [tx_e,tx_n,tx_u];
rx_enu = [rx_e,rx_n,rx_u];

ts = s0_enu-tx_enu; sr = rx_enu-s0_enu;

%convert s1-s4 to the same local ENU
s1 = lla2enu([lat(1) lon(2) ele(1,2)],s0,'ellipsoid');  %north
s2 = lla2enu([lat(3) lon(2) ele(3,2)],s0,'ellipsoid');  %south
s3 = lla2enu([lat(2) lon(3) ele(2,1)],s0,'ellipsoid');  %east
s4 = lla2enu([lat(2) lon(1) ele(2,3)],s0,'ellipsoid');  %west

%local unit North, East and Up vectors
unit_e = (s3-s4)/norm(s3-s4);
unit_n = (s1-s2)/norm(s1-s2);
unit_u = cross(unit_e,unit_n);

p_1e = dot(ts,unit_e); p_1n = dot(ts,unit_n); p_1u = dot(ts,unit_u);
p_2e = dot(sr,unit_e); p_2n = dot(sr,unit_n); p_2u = dot(sr,unit_u);

term1 = p_1e*p_1e+p_1n*p_1n; term2 = p_2e*p_2e+p_2n*p_2n;

theta_i = atand(sqrt(term1)/abs(p_1u)); theta_s = atand(sqrt(term2)/p_2u);
phi_i = atand(p_1n/p_1e); phi_s = atand(p_2n/p_2e);
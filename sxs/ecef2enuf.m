% this function computes the elevation (theta) and azimuth (phi) angle of a point
% in the object's ENU frame (enuf)
% input:
% 1) P: object's ECEF position vector
% 2) S_ecef: ECEF coordinate of the point to be computed
% output: 
% 1) theta_enuf & phi_enuf: elevation and azimuth angles of S in enuf in degree

function [theta_enuf,phi_enuf] = ecef2enuf(P,S_ecef)

%P = [-4593021.50000000	608280.500000000	-4370184.50000000];
%S_ecef = [-4590047.30433596	610685.547457113	-4371634.83935421];

%define wgs84 ellipsoid
wgs84 = wgs84Ellipsoid;

P_lla = ecef2lla(P);

[S_east,S_north,S_up] = ecef2enu(S_ecef(1),S_ecef(2),S_ecef(3),P_lla(1),P_lla(2),P_lla(3),wgs84);
[phi_enuf,theta_enuf1,~] = cart2sph(S_east,S_north,S_up);

phi_enuf = rad2deg(phi_enuf);
theta_enuf1 = rad2deg(theta_enuf1);

theta_enuf = 90-theta_enuf1;
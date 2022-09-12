% this function computes the elevation (theta) and azimuth (phi) angle of a point
% in the object's orbit reference frame (orf)
% Input (all vectors are row vectors):
% 1) P & V: object's ECEF position (P) and velocity (V) vectors
% 2) S_ecef: ECEF coordinate of the point to be computed (S_ecef)
% Output: 
% 1) theta_orf & phi_orf: polar and azimuth angles of S in SV's orf in degree
% 2) S_orf: coordinate of S in orf S_orf

function [theta_orf,phi_orf] = ecef2orf(P,V,S_ecef)

P = P';
V = V';
S_ecef = S_ecef';

u_ecef = S_ecef-P;                      %vector from P to S

theta_e = 7.2921158553e-5;              %earth rotation rate, rad/s
W_e = [0 0 theta_e]';                   %earth rotation vector
Vi = V+cross(W_e,P);                    %SC ECEF inertial velocity vector

%define orbit reference frame - unit vectors
y_orf = cross(-1*P,Vi)/norm(cross(-1*P,Vi));
z_orf = -1*P/norm(P);
x_orf = cross(y_orf,z_orf);

%transformation matrix
T_orf = [x_orf'; y_orf'; z_orf'];

S_orf = T_orf*u_ecef;

%elevation and azimuth angles
theta_orf = acosd(S_orf(3)/(norm(S_orf)));
phi_orf = atan2d(S_orf(2),S_orf(1));

if phi_orf < 0
    phi_orf = 360+phi_orf;
end
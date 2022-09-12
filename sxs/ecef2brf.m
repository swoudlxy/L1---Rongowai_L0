% this function computes the elevation (theta) and azimuth (phi) angle of a
% ecef vector in the objects's body reference frame (brf)
% Input: 
% 1) P, V: object's ecef position vector
% 2) SC_att: object's attitude (Euler angle) in the sequence of
% roll, pitch, yaw, in degrees
% 3) S_ecef: ecef coordinate of the point to be computed
% Output: 
% 1) theta_brf: elevation angle of S in the SC's brf in degree
% 2) phi_brf: azimuth angle of S in the SC's brf in degree

function [theta_brf,phi_brf] = ecef2brf(P,V,S_ecef,SC_att)

P = P'; V = V'; S_ecef = S_ecef';

phi     = deg2rad(SC_att(1));       %roll
theta   = deg2rad(SC_att(2));       %pitch
psi     = deg2rad(SC_att(3));       %yaw

u_ecef = S_ecef-P;                  %vector from P to S

%define heading frame - unit vectors
y_hrf = cross(-1*P,V)/norm(cross(-1*P,V));
z_hrf = -1*P/norm(-1*P);
x_hrf = cross(y_hrf,z_hrf);

T_hrf = [x_hrf';y_hrf';z_hrf'];

%S in hrf
S_hrf = T_hrf*u_ecef;

%construct aircraft's attitude matrix
Rx_phi =    [   1           0               0               ;
                0           cos(phi)        sin(phi)        ;
                0           -1*sin(phi)     cos(phi)]       ;

Ry_theta =  [   cos(theta)  0               -1*sin(theta)   ;
                0           1               0               ;
                sin(theta)  0               cos(theta)]     ;

Rz_psi =    [   cos(psi)    sin(psi)        0               ;
                -1*sin(psi) cos(psi)        0               ;
                0           0               1]              ;

R = Ry_theta*Rx_phi*Rz_psi;         %transformation matrix

S_brf = R*S_hrf;

theta_brf = acosd(S_brf(3)/(norm(S_brf)));
phi_brf = atan2d(S_brf(2),S_brf(1));

if phi_brf < 0
    phi_brf = 360+phi_brf;
end
% this function computes the elevation (theta) and azimuth (phi) angle of a point
% in the SC (spacecraft) body reference frame (brf)
% input (u_orf in column vector): 
% 1) u_orf: orf coordinate of the point to be computed
% 2) att: SC attitude matrix - Eular angles (roll, pitch, yaw) - in radian
% output: 
% 1) polar and azimuth angles of S in the SC's orbit frame in degree
% 2) urf coordinate of the point

function [theta_brf,phi_brf,u_brf] = orf2brf(att,u_orf)

u_orf = u_orf';

%function starts here
phi     = att(1);               %roll
theta   = att(2);               %pitch
psi     = att(3);               %yaw

%construct rotation marix
Rx_phi =    [   1           0               0               ;
                0           cos(phi)        sin(phi)        ;
                0           -1*sin(phi)     cos(phi)]       ;

Ry_theta =  [   cos(theta)  0               -1*sin(theta)   ;
                0           1               0               ;
                sin(theta)  0               cos(theta)]     ;

Rz_psi =    [   cos(psi)    sin(psi)        0               ;
                -1*sin(psi) cos(psi)        0               ;
                0           0               1]              ;

R = Ry_theta*Rx_phi*Rz_psi;     %transformation matrix

u_brf = R*u_orf;

theta_brf = acosd(u_brf(3)/(norm(u_brf)));
phi_brf = atan2d(u_brf(2),u_brf(1));

if phi_brf < 0
    phi_brf = 360+phi_brf;
end
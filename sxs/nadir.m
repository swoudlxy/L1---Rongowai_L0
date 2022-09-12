% This function, based on the WGS84 model, computes the ECEF coordinate
% of the nadir point (n) of a in-space point (m) on a WGS84 model

function [s_ecef]= nadir(m_ecef)

%define WGS84
wgs84 = wgs84Ellipsoid;     %WGS84 in metres
a = wgs84.SemimajorAxis;
invf = wgs84.InverseFlattening;
f = 1/invf;
e2 = 2*f-f*f;               %eccentricity square WGS84

%calculate nadir point
thetaD = asind(m_ecef(3)/norm(m_ecef));     %latitude degrees
cost = cosd(thetaD);
cost2 = cost*cost;

term1 = 1-e2;
term2 = 1-e2*cost2;
r = a*sqrt(term1/term2);                    %lat-dependent Earth radius

s_ecef = r*m_ecef/norm(m_ecef);             %s: nadir of m on WGS84
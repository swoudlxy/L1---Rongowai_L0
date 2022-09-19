% This function determines if the RT vector has intersections
% with the WGS84 ellipsoid (LOS existence)
% input: tx and rx locations in ECEF frame
% output: flag to indicate if LOS exists between tx and rx

function [flag] = los_status(tx,rx)

% define WGS84
wgs84 = wgs84Ellipsoid;

a = wgs84.SemimajorAxis;
b = a;
c = wgs84.SemiminorAxis;

% rx for NGRx, tx for satellite, given in ECEF-XYZ
T_ecef = tx./[a,b,c];   % sat pos vector
R_ecef = rx./[a,b,c];	% ngrx pos vector

RT = T_ecef-R_ecef;     % RT vector
RT_unit = RT/norm(RT);  % unit vector of RT

% determine if LOS exists (flag = 1)
A = norm(RT_unit)*norm(RT_unit);
B = 2*dot(R_ecef,RT_unit); B2 = B*B;
C = norm(R_ecef)*norm(R_ecef)-1;

t1 = B2-4*A*C;
t2 = -B+sqrt(t1)/(2*A);
t3 = -B-sqrt(t1)/(2*A);

if t1 < 0
    flag = 1;
elseif (t2 < 0) && (t3 < 0)
    flag = 1;
else
    flag = 0;
end

end
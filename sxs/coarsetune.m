% this function computes the SP on a pure WGS84 datum based on
% Inputs:
% 1) TX_xyz: ECEF coordinate of the TX
% 2) RX_xyz: ECEF coordinate of the RX
% Outputs:
% 1) SP_xyz, SP_lla: ECEF and LLA coordinate of a SP on a pure WGS84 datum

function [SP_xyz_coarse,SP_lla_coarse] = coarsetune(TX_xyz,RX_xyz)

%find coarse SP using Fibonacci sequence
s1 = ite(TX_xyz,RX_xyz);      
s2 = ecef2lla(s1);

%longitude adjustment
s_lon = s2(2);
if s_lon < 0
    s_lon = s_lon+360;
elseif s_lon > 360
    s_lon = s_lon-360;
end

s2(2) = s_lon;

%save to ecef and lla coordinates
SP_xyz_coarse = s1;
SP_lla_coarse = s2;
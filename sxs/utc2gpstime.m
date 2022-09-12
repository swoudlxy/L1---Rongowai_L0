% this function converts from a utc time to a gps week and sec
% Inputs: UTC time in matrix [year,month,day,hour,minute,sec]
% Outputs: gps week and gps seconds of the week
% Note the utc_offset needs to be adjusted according to here:
% https://en.racelogic.support/LabSat_GNSS_Simulators/Knowledge_Base/General_Info/LabSat_Leap_Second_Guide

function [gps_week,gps_sec] = utc2gpstime(year,month,day,hour,min,sec)
    
% Number of days into the year at the start of each month (ignoring leap years).
days_in_year = [0,31,59,90,120,151,181,212,243,273,304,334,365];

% year referenced from 1980
ye = year - 1980;

% compute num of leap days
leap_days = ye/4 + 1;
if (mod(ye,4) == 0) && (month <= 2)
    leap_days = leap_days - 1;
end

% days elapsed since midnight jan 5, 1980.
de = floor(ye*365 + days_in_year(month) + day + leap_days - 6);

% leap Seconds, good after 1998
% GPS time is ahead of UTC by this many Leap Seconds
utc_offset = 18;        % 2016/12/31 - XXX

if year < 1999
    sprintf('Leap seconds invalid for this year!\n');
    gps_week = 0;
    gps_sec = 0;
    return
end

% Convert time to GPS weeks and seconds.
gps_week = floor(de/7.0);
gps_sec = mod(de,7)*86400.0 + hour*3600.0 + min*60.0 + sec + utc_offset;

% Adjust GPS weeks/seconds to guarantee that seconds is in the range 0-604800.0. */
if gps_sec < 0.0
    gps_week = gps_week-1;
    gps_sec = gps_sec+64800.0;
end

if gps_sec >= 604800.0        
    gps_week = gps_week+1;
    gps_sec = gps_sec-64800.0;
end

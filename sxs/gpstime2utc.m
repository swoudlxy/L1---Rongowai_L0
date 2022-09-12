% this function converts from a gps time to a utc time
% Inputs: gps week and gps seconds of the week
% Outputs: UTC time in matrix [year,month,day,hour,minute,sec]
% Note the utc_offset needs to be adjusted according to here:
% https://en.racelogic.support/LabSat_GNSS_Simulators/Knowledge_Base/General_Info/LabSat_Leap_Second_Guide

function utc_time = gpstime2utc(gpsweek,gpssecs)

days_in_year = [0,31,59,90,120,151,181,212,243,273,304,334,365];
days_in_leap_year = [0,31,60,91,121,152,182,213,244,274,305,335,366];
    
% leap days since 1980
leap_day_table = [1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0, ...
    1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0];

% adjust for utc offset
utc_offset = 18;        % 2016/12/31 - XXX
secs = gpssecs - utc_offset;

% time as integer and fractional seconds since midnight Jan 5 1980
int_secs = floor(secs);
fract_secs = secs - int_secs;
secs_1980 = 604800*gpsweek + int_secs;

% express current GPS time as elapsed UTC days and seconds.
day_number = floor(secs_1980/86400);
sec_of_day = secs_1980 - day_number*86400;

leap_days1 = floor((day_number+1406) / 1461);

% calculate UTC year
total_years = floor((day_number + 5 - leap_days1) / 365);
year = 1980 + total_years;

% day of utc year
leap_days2 = sum(leap_day_table(1:total_years+1));
day_of_utc_year = (day_number + 5 - 365*total_years - leap_days2);

% determine month and day
month = -1;

if rem(year,4) ~= 0    
    % not a leap year
    for i = 1:13
        if day_of_utc_year > days_in_year(i)
            month = i;
            if month == 13
                month = 12;
            end

            day = day_of_utc_year - days_in_year(i) + 1;
        end
    end
    
else
    % a leap year
    for i = 1:13
        if day_of_utc_year > days_in_leap_year(i)
            month = i;
            if month == 13
                month = 12;
            end

            day = day_of_utc_year - days_in_leap_year(i) + 1;
        end
    end
    
end 

% calculate utc hour
hour = floor(sec_of_day/3600);

if hour > 23
    hour = 23;
end

% calculate utc minute
minute = floor((sec_of_day - hour*3600) / 60);

if minute > 59
    minute = 59;
end

% calculate utc seconds
sec = (sec_of_day - hour*3600 - minute*60) + fract_secs;

utc_time = [year,month,day,hour,minute,sec];

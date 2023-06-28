% this function retrieves GNSS antenna pattern at GPS/GAL L1/E1 band
% this function outputs the gain LUT for a specific polarisation manner

function [azim_deg,elev_deg,gain_data] = get_ant_pattern(ant_filename)

fid = fopen(ant_filename,'r');

% update precision type 27 June
azim_min = fread(fid,1,'uint16');
azim_max = fread(fid,1,'uint16');

elev_max = fread(fid,1,'uint16');
elev_min = fread(fid,1,'uint16');

res = fread(fid,1,'double');

ant_data = fread(fid,3601*1201,'double');

fclose(fid);

azim_deg = azim_min:res:azim_max;
elev_deg = elev_max:-1*res:elev_min;
gain_data = reshape(ant_data,[3601,1201]);
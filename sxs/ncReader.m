% this is a simple netCDF reader
% this reader first reads all the variable names from netCDF
% and then retrieve the variables under each name and save to
% the workspace
clear
clc

% define file and retrieve netCDF variables, attributes and comments
filename = '../out/L1_sample1.nc';
L1_sample = ncinfo(filename);

% retrieve variables and save to workspace
var_name = L1_sample.Variables;
L = length(var_name);

for l = 1:L

    var_name1 = var_name(l).Name;

    a = double(ncread(filename,var_name1));

    eval([var_name1 '= a;']);

end

clear a var_name var_name1 L;                   % clear variables do not belong to netCDF
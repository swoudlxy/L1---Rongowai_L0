% load data
clear
clc

load('../out/L1_postCal.mat')
L1_dict = readtable('../dat/L1_dict_final.xlsx');
L1_dict = string(table2cell(L1_dict));

field_names = fieldnames(L1_postCal);
L = length(field_names);

netCDF_name = '../out/L1_sample1.nc';

%% Part 1
clc

ncid = netcdf.create(netCDF_name,'NETCDF4');

for l = 1:L

    field_name1 = string(field_names(l));
    index1 = strcmp(L1_dict(:,1),field_name1);

    netCDF_field_name = L1_dict(index1,1);
    unit = L1_dict(index1,4);
    dimension = L1_dict(index1,5);

    % global attributes
    if strcmp(unit,'<none>') && strcmp(dimension,'<none>')
        att_value = getfield(L1_postCal,field_name1);
        netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),netCDF_field_name,string(att_value));
    end

end

netcdf.close(ncid);

%% Part 2
clc

% retrieve dimensions for netCDF
sample = length(L1_postCal.ddm_timestamp_utc);
ddm = 20;
delay = 40;
doppler = 5;

for l = 1:30

    field_name1 = string(field_names(l));
    val_value = getfield(L1_postCal,field_name1);

    index1 = strcmp(L1_dict(:,1),field_name1);
    netCDF_field_name = L1_dict(index1,1);
    long_name = L1_dict(index1,2);
    data_type = L1_dict(index1,3);
    unit = L1_dict(index1,4);
    dimension = L1_dict(index1,5);
    comment = L1_dict(index1,6);

    if ~strcmp(unit,'<none>') && strcmp(dimension,'<none>')
        nccreate(netCDF_name,netCDF_field_name);
        ncwrite(netCDF_name,netCDF_field_name,val_value);

        ncwriteatt(netCDF_name,netCDF_field_name,'long name',long_name);
        ncwriteatt(netCDF_name,netCDF_field_name,'units',unit);
        ncwriteatt(netCDF_name,netCDF_field_name,'comment',comment);
    end

end

for l = 31:L

field_name1 = string(field_names(l));
val_value = getfield(L1_postCal,field_name1);

index1 = strcmp(L1_dict(:,1),field_name1);
netCDF_field_name = L1_dict(index1,1);
long_name = L1_dict(index1,2);
data_type = L1_dict(index1,3);
unit = L1_dict(index1,4);
dimension = L1_dict(index1,5);
comment = L1_dict(index1,6);

    if strcmp(dimension,'sample')
        val_value = cast(val_value,data_type);

        nccreate(netCDF_name,netCDF_field_name, ...
                 'Dimensions',{'sample',sample}, ...
                 'FillValue','disable');
        ncwrite(netCDF_name,netCDF_field_name,val_value);

        ncwriteatt(netCDF_name,netCDF_field_name,'long name',long_name);
        ncwriteatt(netCDF_name,netCDF_field_name,'units',unit);
        ncwriteatt(netCDF_name,netCDF_field_name,'comment',comment);

    elseif strcmp(dimension,'ddm')
        val_value = cast(val_value,data_type);

        nccreate(netCDF_name,netCDF_field_name, ...
                 'Dimensions',{'ddm',ddm}, ...
                 'FillValue','disable');
        ncwrite(netCDF_name,netCDF_field_name,val_value);

        ncwriteatt(netCDF_name,netCDF_field_name,'long name',long_name);
        ncwriteatt(netCDF_name,netCDF_field_name,'units',unit);
        ncwriteatt(netCDF_name,netCDF_field_name,'comment',comment);

    elseif strcmp(dimension,'sample, ddm')
        val_value = cast(val_value,data_type);

        nccreate(netCDF_name,netCDF_field_name, ...
                 'Dimensions',{'ddm',ddm,'sample',sample}, ...
                 'FillValue','disable');
        ncwrite(netCDF_name,netCDF_field_name,val_value);

        ncwriteatt(netCDF_name,netCDF_field_name,'long name',long_name);
        ncwriteatt(netCDF_name,netCDF_field_name,'units',unit);
        ncwriteatt(netCDF_name,netCDF_field_name,'comment',comment);

    elseif strcmp(dimension,'sample, ddm, delay, doppler')
        val_value = cast(val_value,data_type);

        nccreate(netCDF_name,netCDF_field_name, ...
                 'Dimensions',{'doppler',doppler,'delay',delay,'ddm',ddm,'sample',sample}, ...
                 'FillValue','disable');
        ncwrite(netCDF_name,netCDF_field_name,val_value);

        ncwriteatt(netCDF_name,netCDF_field_name,'long name',long_name);
        ncwriteatt(netCDF_name,netCDF_field_name,'units',unit);
        ncwriteatt(netCDF_name,netCDF_field_name,'comment',comment);

    elseif strcmp(dimension,'sample, ddm, delay')
        val_value = cast(val_value,data_type);
        nccreate(netCDF_name,netCDF_field_name, ...
                 'Dimensions',{'delay',delay,'ddm',ddm,'sample',sample}, ...
                 'FillValue','disable');
        ncwrite(netCDF_name,netCDF_field_name,val_value);

        ncwriteatt(netCDF_name,netCDF_field_name,'long name',long_name);
        ncwriteatt(netCDF_name,netCDF_field_name,'units',unit);
        ncwriteatt(netCDF_name,netCDF_field_name,'comment',comment);
    end

end

%{
for l = 20:L

    field_name1 = string(field_names(l));
    val_value = getfield(L1_postCal,field_name1);

    index1 = strcmp(L1_dict(:,1),field_name1);
    netCDF_field_name = L1_dict(index1,1);
    long_name = L1_dict(index1,2);
    data_type = L1_dict(index1,3);
    unit = L1_dict(index1,4);
    dimension = L1_dict(index1,5);
    comment = L1_dict(index1,6);

    if ~strcmp(unit,'<none>') && strcmp(dimension,'<none>')
        nccreate(netCDF_name,netCDF_field_name);
        ncwrite(netCDF_name,netCDF_field_name,val_value);

        ncwriteatt(netCDF_name,netCDF_field_name,'long name',long_name);
        ncwriteatt(netCDF_name,netCDF_field_name,'units',unit);
        ncwriteatt(netCDF_name,netCDF_field_name,'comment',comment);
 
    elseif strcmp(dimension,'sample')
        %val_value = cast(val_value,data_type);

        nccreate(netCDF_name,netCDF_field_name, ...
                 'Dimensions',{'sample',sample}, ...
                 'FillValue','disable');
        ncwrite(netCDF_name,netCDF_field_name,val_value);

        ncwriteatt(netCDF_name,netCDF_field_name,'long name',long_name);
        ncwriteatt(netCDF_name,netCDF_field_name,'units',unit);
        ncwriteatt(netCDF_name,netCDF_field_name,'comment',comment);

    elseif strcmp(dimension,'sample, ddm')
        %val_value = cast(val_value,data_type);

        nccreate(netCDF_name,netCDF_field_name, ...
                 'Dimensions',{'ddm',ddm,'sample',sample}, ...
                 'FillValue','disable');
        ncwrite(netCDF_name,netCDF_field_name,val_value);

        ncwriteatt(netCDF_name,netCDF_field_name,'long name',long_name);
        ncwriteatt(netCDF_name,netCDF_field_name,'units',unit);
        ncwriteatt(netCDF_name,netCDF_field_name,'comment',comment);

    elseif strcmp(dimension,'sample, ddm, delay, doppler')
        %val_value = cast(val_value,data_type);

        nccreate(netCDF_name,netCDF_field_name, ...
                 'Dimensions',{'doppler',doppler,'delay',delay,'ddm',ddm,'sample',sample}, ...
                 'FillValue','disable');
        ncwrite(netCDF_name,netCDF_field_name,val_value);

        ncwriteatt(netCDF_name,netCDF_field_name,'long name',long_name);
        ncwriteatt(netCDF_name,netCDF_field_name,'units',unit);
        ncwriteatt(netCDF_name,netCDF_field_name,'comment',comment);

    elseif strcmp(dimension,'sample, ddm, delay')
        %val_value = cast(val_value,data_type);
        nccreate(netCDF_name,netCDF_field_name, ...
                 'Dimensions',{'delay',delay,'ddm',ddm,'sample',sample}, ...
                 'FillValue','disable');
        ncwrite(netCDF_name,netCDF_field_name,val_value);

        ncwriteatt(netCDF_name,netCDF_field_name,'long name',long_name);
        ncwriteatt(netCDF_name,netCDF_field_name,'units',unit);
        ncwriteatt(netCDF_name,netCDF_field_name,'comment',comment);

    end
end
%}

%%
%netCDF_field_name = 'sample';
%val_value = L1_postCal.sample;
%{
nccreate(netCDF_name,netCDF_field_name, ...
    'Dimensions',{'sample',sample}, ...
    'FillValue','disable');
ncwrite(netCDF_name,netCDF_field_name,val_value);

ncwriteatt(netCDF_name,netCDF_field_name,'long name',long_name);
ncwriteatt(netCDF_name,netCDF_field_name,'units',unit);
ncwriteatt(netCDF_name,netCDF_field_name,'comment',comment);
%}
L1_sample = ncinfo('../out/L1_sample1.nc');
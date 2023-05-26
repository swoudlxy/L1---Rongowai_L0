function sample_info = get_netcdf(L1_netCDF_name,L1_dict_name,L1_postCal)

L1_dict = readtable(L1_dict_name);
L1_dict = string(table2cell(L1_dict));

L1_postCal = rmfield(L1_postCal,'A_eff_all');
field_names = fieldnames(L1_postCal);
%field_names(117) = [];                  % debug only

L = length(field_names);

ncid = netcdf.create(L1_netCDF_name,'NETCDF4');

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

% retrieve dimensions for netCDF
num_sample = length(L1_postCal.ddm_timestamp_utc);
num_ddm = 20;
delay = 40;
doppler = 5;

for l = 1:30

    field_name1 = string(field_names(l));
    val_value = getfield(L1_postCal,field_name1);

    index1 = strcmp(L1_dict(:,1),field_name1);
    netCDF_field_name = L1_dict(index1,1);
    long_name = L1_dict(index1,2);
    unit = L1_dict(index1,4);
    dimension = L1_dict(index1,5);
    comment = L1_dict(index1,6);

    if ~strcmp(unit,'<none>') && strcmp(dimension,'<none>')
        nccreate(L1_netCDF_name,netCDF_field_name);
        ncwrite(L1_netCDF_name,netCDF_field_name,val_value);

        ncwriteatt(L1_netCDF_name,netCDF_field_name,'long name',long_name);
        ncwriteatt(L1_netCDF_name,netCDF_field_name,'units',unit);
        ncwriteatt(L1_netCDF_name,netCDF_field_name,'comment',comment);

    elseif strcmp(dimension,'sample')
        nccreate(L1_netCDF_name,netCDF_field_name, ...
                 'Dimensions',{'sample',num_sample}, ...
                 'FillValue','disable');
        ncwrite(L1_netCDF_name,netCDF_field_name,val_value);

        ncwriteatt(L1_netCDF_name,netCDF_field_name,'long name',long_name);
        ncwriteatt(L1_netCDF_name,netCDF_field_name,'units',unit);
        ncwriteatt(L1_netCDF_name,netCDF_field_name,'comment',comment);
    
    end

end

sample_value = (0:1:num_sample-1)';
nccreate(L1_netCDF_name,'sample_index', ...
         'Dimensions',{'sample',num_sample}, ...
         'FillValue','disable');
ncwrite(L1_netCDF_name,'sample_index',sample_value);

for l = 32:L

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
        nccreate(L1_netCDF_name,netCDF_field_name, ...
                 'Dimensions',{'sample',num_sample}, ...
                 'FillValue','disable');
        ncwrite(L1_netCDF_name,netCDF_field_name,val_value);

        ncwriteatt(L1_netCDF_name,netCDF_field_name,'long name',long_name);
        ncwriteatt(L1_netCDF_name,netCDF_field_name,'units',unit);
        ncwriteatt(L1_netCDF_name,netCDF_field_name,'comment',comment);

    elseif strcmp(dimension,'ddm')
        val_value = cast(val_value,data_type);

        nccreate(L1_netCDF_name,netCDF_field_name, ...
                 'Dimensions',{'ddm',num_ddm}, ...
                 'FillValue','disable');
        ncwrite(L1_netCDF_name,netCDF_field_name,val_value);

        ncwriteatt(L1_netCDF_name,netCDF_field_name,'long name',long_name);
        ncwriteatt(L1_netCDF_name,netCDF_field_name,'units',unit);
        ncwriteatt(L1_netCDF_name,netCDF_field_name,'comment',comment);

    elseif strcmp(dimension,'sample, ddm')
        val_value = cast(val_value,data_type);

        nccreate(L1_netCDF_name,netCDF_field_name, ...
                 'Dimensions',{'ddm',num_ddm,'sample',num_sample}, ...
                 'FillValue','disable');
        ncwrite(L1_netCDF_name,netCDF_field_name,val_value);

        ncwriteatt(L1_netCDF_name,netCDF_field_name,'long name',long_name);
        ncwriteatt(L1_netCDF_name,netCDF_field_name,'units',unit);
        ncwriteatt(L1_netCDF_name,netCDF_field_name,'comment',comment);

    elseif strcmp(dimension,'sample, ddm, delay, doppler')
        val_value = cast(val_value,data_type);

        nccreate(L1_netCDF_name,netCDF_field_name, ...
                 'Dimensions',{'doppler',doppler,'delay',delay,'ddm',num_ddm,'sample',num_sample}, ...
                 'FillValue','disable');
        ncwrite(L1_netCDF_name,netCDF_field_name,val_value);

        ncwriteatt(L1_netCDF_name,netCDF_field_name,'long name',long_name);
        ncwriteatt(L1_netCDF_name,netCDF_field_name,'units',unit);
        ncwriteatt(L1_netCDF_name,netCDF_field_name,'comment',comment);

    elseif strcmp(dimension,'sample, ddm, delay')
        val_value = cast(val_value,data_type);
        nccreate(L1_netCDF_name,netCDF_field_name, ...
                 'Dimensions',{'x',1,'delay',delay,'ddm',num_ddm,'sample',num_sample}, ...
                 'FillValue','disable');
        ncwrite(L1_netCDF_name,netCDF_field_name,val_value);

        ncwriteatt(L1_netCDF_name,netCDF_field_name,'long name',long_name);
        ncwriteatt(L1_netCDF_name,netCDF_field_name,'units',unit);
        ncwriteatt(L1_netCDF_name,netCDF_field_name,'comment',comment);

    end

end

sample_info = ncinfo(L1_netCDF_name);
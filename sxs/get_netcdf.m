function sample_info = get_netcdf(L1_netCDF_name,L1_dict_name,L1_postCal)

L1_dict = readtable(L1_dict_name);
L1_dict = string(table2cell(L1_dict));

field_names = fieldnames(L1_postCal);

L = length(field_names);

ncid = netcdf.create(L1_netCDF_name,'NETCDF4');

for l = 1:42

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

for l = 1:49

    field_name1 = string(field_names(l));
    val_value = getfield(L1_postCal,field_name1);

    index1 = strcmp(L1_dict(:,1),field_name1);
    netCDF_field_name = L1_dict(index1,1);
    long_name = L1_dict(index1,2);
    dtype = L1_dict(index1,3);
    unit = L1_dict(index1,4);
    dimension = L1_dict(index1,5);
    comment = L1_dict(index1,6);

    if ~strcmp(unit,'<none>') && strcmp(dimension,'<none>')
        nccreate(L1_netCDF_name,netCDF_field_name,"Datatype",dtype);
        ncwrite(L1_netCDF_name,netCDF_field_name,val_value);

        ncwriteatt(L1_netCDF_name,netCDF_field_name,'long_name',long_name);
        ncwriteatt(L1_netCDF_name,netCDF_field_name,'units',unit);
        ncwriteatt(L1_netCDF_name,netCDF_field_name,'comment',comment);

    elseif strcmp(dimension,'sample')
        nccreate(L1_netCDF_name,netCDF_field_name, ...
                 'Dimensions',{'sample',num_sample}, ...
                 'FillValue','disable','Datatype',dtype);
        ncwrite(L1_netCDF_name,netCDF_field_name,val_value);

        ncwriteatt(L1_netCDF_name,netCDF_field_name,'long_name',long_name);

        % compliance check updat - 30 June
        if strcmp(field_name1,'ddm_timestamp_utc') || strcmp(field_name1,'pvt_timestamp_utc')
            ncwriteatt(L1_netCDF_name,netCDF_field_name,'standard_name','time');
            ncwriteatt(L1_netCDF_name,netCDF_field_name,'calendar','gregorian');

            start_time = L1_postCal.time_coverage_start;

            ncwriteatt(L1_netCDF_name,netCDF_field_name,'units', ...
                join(['seconds since' string(datetime(start_time,'format','yyyy-MM-dd HH:mm:ss'))]));

        else
            ncwriteatt(L1_netCDF_name,netCDF_field_name,'units',unit);

        end

        ncwriteatt(L1_netCDF_name,netCDF_field_name,'comment',comment);
    
    end

end

% field 50
% compliance check update 30 June
sample_value = (0:1:num_sample-1)';
nccreate(L1_netCDF_name,'sample_index', ...
         'Dimensions',{'sample',num_sample}, ...
         'FillValue','disable','Datatype',dtype);
ncwrite(L1_netCDF_name,'sample_index',sample_value);

ncwriteatt(L1_netCDF_name,'sample_index','long_name','Sample index');
ncwriteatt(L1_netCDF_name,'sample_index','units','1');
ncwriteatt(L1_netCDF_name,'sample_index','comment', ...
    join(['The netCDF coordinate variable associated with the sample dimension, ' ...
    'which enumerates the zero-justified index range of the DDM time instants contained in the file.']));

for l = 51:L

    field_name1 = string(field_names(l));
    val_value = getfield(L1_postCal,field_name1);

    index1 = strcmp(L1_dict(:,1),field_name1);
    netCDF_field_name = L1_dict(index1,1);

    long_name = L1_dict(index1,2);
    dtype = L1_dict(index1,3);
    unit = L1_dict(index1,4);
    dimension = L1_dict(index1,5);
    comment = L1_dict(index1,6);

    if strcmp(dimension,'sample')
        nccreate(L1_netCDF_name,netCDF_field_name, ...
                 'Dimensions',{'sample',num_sample}, ...
                 'FillValue','disable','Datatype',dtype);
        ncwrite(L1_netCDF_name,netCDF_field_name,val_value);

        ncwriteatt(L1_netCDF_name,netCDF_field_name,'long_name',long_name);
        ncwriteatt(L1_netCDF_name,netCDF_field_name,'units',unit);

        % compliance check update - 30 June
        if strcmp(field_name1,'ac_lat')
            ncwriteatt(L1_netCDF_name,netCDF_field_name,'standard_name','latitude');

        elseif strcmp(field_name1,'ac_lon')
            ncwriteatt(L1_netCDF_name,netCDF_field_name,'standard_name','longitude');

        end

        if ~strcmp(unit,'<none>')
            ncwriteatt(L1_netCDF_name,netCDF_field_name,'units',unit);
        end
        
        ncwriteatt(L1_netCDF_name,netCDF_field_name,'comment',comment);

    elseif strcmp(dimension,'ddm')
        val_value = cast(val_value,dtype);

        nccreate(L1_netCDF_name,netCDF_field_name, ...
                 'Dimensions',{'ddm',num_ddm}, ...
                 'FillValue','disable','Datatype',dtype);
        ncwrite(L1_netCDF_name,netCDF_field_name,val_value);

        ncwriteatt(L1_netCDF_name,netCDF_field_name,'long_name',long_name);
        ncwriteatt(L1_netCDF_name,netCDF_field_name,'units',unit);
        ncwriteatt(L1_netCDF_name,netCDF_field_name,'comment',comment);

    elseif strcmp(dimension,'sample, ddm')
        val_value = cast(val_value,dtype);

        nccreate(L1_netCDF_name,netCDF_field_name, ...
                 'Dimensions',{'ddm',num_ddm,'sample',num_sample}, ...
                 'FillValue','disable','Datatype',dtype);
        ncwrite(L1_netCDF_name,netCDF_field_name,val_value);

        ncwriteatt(L1_netCDF_name,netCDF_field_name,'long_name',long_name);

        % compliance check update - 30 June
        if strcmp(field_name1,'sp_lat')
            ncwriteatt(L1_netCDF_name,netCDF_field_name,'standard_name','latitude');

        elseif strcmp(field_name1,'sp_lon')
            ncwriteatt(L1_netCDF_name,netCDF_field_name,'standard_name','longitude');

        end

        if ~strcmp(unit,'<none>')
            ncwriteatt(L1_netCDF_name,netCDF_field_name,'units',unit);
        end
        
        ncwriteatt(L1_netCDF_name,netCDF_field_name,'comment',comment);
        
    elseif strcmp(dimension,'sample, ddm, delay, doppler')
        val_value = cast(val_value,dtype);

        nccreate(L1_netCDF_name,netCDF_field_name, ...
                 'Dimensions',{'doppler',doppler,'delay',delay,'ddm',num_ddm,'sample',num_sample}, ...
                 'FillValue','disable','Datatype',dtype);
        ncwrite(L1_netCDF_name,netCDF_field_name,val_value);

        ncwriteatt(L1_netCDF_name,netCDF_field_name,'long_name',long_name);
        ncwriteatt(L1_netCDF_name,netCDF_field_name,'units',unit);
        ncwriteatt(L1_netCDF_name,netCDF_field_name,'comment',comment);

    elseif strcmp(dimension,'sample, ddm, delay')
        val_value = cast(val_value,dtype);
        nccreate(L1_netCDF_name,netCDF_field_name, ...
                 'Dimensions',{'x',1,'delay',delay,'ddm',num_ddm,'sample',num_sample}, ...
                 'FillValue','disable','Datatype',dtype);
        ncwrite(L1_netCDF_name,netCDF_field_name,val_value);

        ncwriteatt(L1_netCDF_name,netCDF_field_name,'long_name',long_name);
        ncwriteatt(L1_netCDF_name,netCDF_field_name,'units',unit);
        ncwriteatt(L1_netCDF_name,netCDF_field_name,'comment',comment);

    end

end

sample_info = ncinfo(L1_netCDF_name);
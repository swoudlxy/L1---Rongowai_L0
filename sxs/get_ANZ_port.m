% this function maps the rf source port to the ANZ RF port to indicate the
% source of DDM

function ANZ_port = get_ANZ_port(rf_source)

if rf_source == 0
    ANZ_port = 1;       % zenith

elseif rf_source == 4
    ANZ_port = 2;       % nadir LHCP

elseif rf_source == 8
    ANZ_port = 3;       % nadir RHCP

end

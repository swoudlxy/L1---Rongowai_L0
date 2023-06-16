% this function constructs a A_phy LUT as a structure

function A_phy_LUT_all = get_A_phy_LUT(A_phy_LUT_path)
%clear
%clc

%A_phy_LUT_path = '../out/A_phy_LUT_new/';

A_phy_LUT_files = dir([A_phy_LUT_path '*.dat']);
L = length(A_phy_LUT_files);

A_phy_LUT_all = struct;

for l = 1:L

    filename = A_phy_LUT_files(l).name;
    LUT_file = [A_phy_LUT_path filename];

    fid = fopen(LUT_file);

    delay_shift = fread(fid,1,'int8');
    doppler_shift = fread(fid,1,'int8');

    A_phy_LUT = zeros(131,21,13)+nan;

    for j = 1:13

        num_samples = 131*21;

        data = fread(fid,num_samples,'uint32');
        data = reshape(data,[131,21]);
        A_phy_LUT(:,:,j) = data;

    end

    fclose(fid);

    A_phy_LUT_all(l).delay_shift = delay_shift;
    A_phy_LUT_all(l).doppler_shift = doppler_shift;
    A_phy_LUT_all(l).A_phy_LUT = A_phy_LUT;

end
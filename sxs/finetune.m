% This code fine-tunes the coordinate of the initial SP based on the DTU10
% datum thorugh a number of iterative steps.

function [res,min_delay,sx_temp] = finetune(tx_xyz,rx_xyz,sx_lla,L,model)

%main function starts here
%define chip length
c = 299792458;                  %light speed
chip_rate = 1.023e6;            %L1 GPS cps
l_chip = c/chip_rate;           %chip length

num_grid = 11;                  %constant grid matrix 11*11

%find the pixel location
sx_lla_temp = sx_lla;
sx_lat_temp = sx_lla_temp(1);
sx_lon_temp = sx_lla_temp(2);

min_lat = sx_lat_temp-L/2; max_lat = sx_lat_temp+L/2;
min_lon = sx_lon_temp-L/2; max_lon = sx_lon_temp+L/2;

lat_bin = linspace(min_lat,max_lat,num_grid);
lon_bin = linspace(min_lon,max_lon,num_grid);

ele = zeros(num_grid);
delay_chip = zeros(num_grid);

for m = 1:num_grid
    for n = 1:num_grid

        p_lat = lat_bin(m);
        p_lon = lon_bin(n);

        p_ele = mapdtu(p_lat,p_lon,model);

        p_xyz = lla2ecef([p_lat p_lon p_ele]);
        p_delay = pdis(tx_xyz,rx_xyz,p_xyz);
        
        ele(m,n) = p_ele;
        delay_chip(m,n) = p_delay/l_chip;

    end
end

%index of the pixel with minimal reflection path
min_delay = min(min(delay_chip));
[m_i,n_i] = find(delay_chip == min_delay);
m_i = m_i(1);n_i = n_i(1);
sx_temp = [lat_bin(m_i) lon_bin(n_i) ele(m_i,n_i)];

%compute resolution
res = m_lldist([lon_bin((num_grid-1)/2) lon_bin((num_grid-1)/2+1)],[lat_bin((num_grid-1)/2),lat_bin((num_grid-1)/2+1)]);
res = res*1000;                 %resolution in metres

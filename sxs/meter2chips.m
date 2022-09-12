%this function converts from meters to chips
%input: x - distance in meters
%output: y - distance in chips

function y = meter2chips(x)

%define constants
c = 299792458;                      %light speed metre per second
chip_rate = 1.023e6;                %L1 GPS chip-per-second, code modulation frequency
tau_c = 1/chip_rate;                %C/A code chiping period 
l_chip = c*tau_c;                   %chip length

y = x/l_chip;

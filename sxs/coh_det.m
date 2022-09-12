% this function computes the coherency of an input raw-count ddm
% Inputs
% 1)raw ddm measured in counts
% 2)SNR measured in decibels
% Outputs
% 1)coherency ratio (CR)
% 2)coherency state (CS)

function [CR,CS] = coh_det(raw_counts,snr_db)

raw_counts = raw_counts';

peak_power = max(max(raw_counts));
[delay_peak,dopp_peak] = find(raw_counts == peak_power,1);  %row and column of peak power

% thermal noise exclusion
% TODO: the threshold may need to be redefined
thre_coeff = 1.055*exp(-0.193*snr_db);
thre = thre_coeff*peak_power;                               %noise exclusion threshold

index = raw_counts < thre;
raw_counts(index) = 0;
    
%deterimine DDMA range
delay_range = delay_peak-1:delay_peak+1;  delay_min = min(delay_range);   delay_max = max(delay_range);
dopp_range = dopp_peak-2:dopp_peak+2;     dopp_min = min(dopp_range);     dopp_max = max(dopp_range);

%determine if DDMA is within DDM, refine if needed
if delay_min < 2
    delay_range = 1:3;

elseif delay_max > 16
    delay_range = 15:17;

end

if dopp_min < 3
    dopp_range = 1:5;

elseif dopp_max > 9
    dopp_range = 7:11;

end

C_in = sum(raw_counts(delay_range,dopp_range),'all');   %summation of DDMA
C_out = sum(raw_counts,'all')-C_in;                     %summation of DDM excluding DDMA

CR = C_in/C_out;                                        %coherency ratio

if CR >= 2
    CS = 1;

elseif CR < 2
    CS = 0;

end
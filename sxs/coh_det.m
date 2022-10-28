% this function computes the coherency of an input raw-count ddm
% Inputs
% 1)raw ddm measured in counts
% 2)SNR measured in decibels
% Outputs
% 1)coherency ratio (CR)
% 2)coherency state (CS)

function [CR,CS] = coh_det(raw_counts,snr_db)

peak_counts = max(max(raw_counts));
[dopp_peak,delay_peak] = find(raw_counts == peak_counts,1);  % row and column of peak power

% thermal noise exclusion
% TODO: the threshold may need to be redefined
thre_coeff = 1.055*exp(-0.193*snr_db);
thre = thre_coeff*peak_counts;                               % noise exclusion threshold

index = raw_counts < thre;
raw_counts(index) = 0;
    
% deterimine DDMA range
delay_range = delay_peak-1:delay_peak+1;  delay_min = min(delay_range);   delay_max = max(delay_range);
dopp_range = dopp_peak-1:dopp_peak+1;     dopp_min = min(dopp_range);     dopp_max = max(dopp_range);

% determine if DDMA is within DDM, refine if needed
if delay_min < 2
    delay_range = 1:3;

elseif delay_max > 39
    delay_range = 38:40;

end

if dopp_min < 2
    dopp_range = 1:3;

elseif dopp_max > 4
    dopp_range = 3:5;

end

C_in = sum(raw_counts(dopp_range,delay_range),'all');   % summation of DDMA
C_out = sum(raw_counts,'all')-C_in;                     % summation of DDM excluding DDMA

CR = C_in/C_out;                                        % coherency ratio

if CR >= 2
    CS = 1;

else    %elseif CR < 2
    CS = 0;
    
end
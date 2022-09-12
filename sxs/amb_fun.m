% this function computes the ambiguity function
% inputs
% 1) tau_s: delay in seconds
% 2) freq_Hz: Doppler in Hz
% 3) tau_c: chipping period in second, 1/chip_rate
% 4) Ti: coherent integration time in seconds
% output
% 1) chi: ambiguity function, product of Lambda and S

function chi = amb_fun(dtau_s,dfreq_Hz,tau_c,Ti)

det = tau_c*(1+tau_c/Ti);               %discriminant for computing Lambda

%compute Lambda - delay
if abs(dtau_s) <= det
    Lambda = 1-abs(dtau_s)/tau_c;

elseif abs(dtau_s) > det
    Lambda = -tau_c/Ti;

end

%compute S - Doppler
S1 = pi*dfreq_Hz*Ti;

if S1 == 0
    S = 1;

elseif S1 ~= 0
    term1 = sin(S1)/S1;
    term2 = exp(-1i*S1);

    S = term1*term2;
end

%compute complex chi
chi = Lambda*S;
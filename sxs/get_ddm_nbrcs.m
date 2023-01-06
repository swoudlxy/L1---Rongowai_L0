% this function computes the normalised BRCS, leading edge slope and
% trailing edge slope of a DDM
% Inputs
% 1) brcs: bistatic radar cross section
% 2) A_eff: effective scattering area
% 3) sp_delay/doppler_bin_float: floating specular point bin
% Outputs:
% 1) nbrcs - structure: NBRCS and NBRCS effective scattering area 
% 3) LES - structure: LES slope and LES effective scattering area
% 4) TES - structure: TES slope

function [nbrcs,LES,TES] = get_ddm_nbrcs(brcs,A_eff,sx)

%brcs = brcs1; A_eff = A_eff1; sx = sx1;

sx_delay_bin_float = sx.sx_delay_bin+1;
sx_doppler_bin_float = sx.sx_doppler_bin+1;

% define DDMA region size
% TO DO: this region size needs to be refined after getting real data
num_delay_bins_ddma = 3;
num_doppler_bins_ddma = 3;

% perform only when sp bins are within DDM region
if (sx_delay_bin_float < 37) && (sx_delay_bin_float > 2) && ...
    (sx_doppler_bin_float < 4) && (sx_doppler_bin_float > 2)
    
    t0 = floor(sx_delay_bin_float);
    dt = sx_delay_bin_float-t0;

    f0 = floor(sx_doppler_bin_float);

    nbrcs_delay_range = t0:t0+num_delay_bins_ddma-1;
    nbrcs_doppler_range = f0-floor(num_doppler_bins_ddma/2):f0+floor(num_doppler_bins_ddma/2);

    coeff1 = 1-dt;
    term1 = [brcs(f0-1,t0)      brcs(f0,t0)     brcs(f0+1,t0)];

    coeff2 = dt;
    term2 = [brcs(f0-1,t0+3)    brcs(f0,t0+3)   brcs(f0+1,t0+3)];

    term3 = [brcs(f0-1,t0+1)    brcs(f0,t0+1)   brcs(f0+1,t0+1), ...
             brcs(f0-1,t0+2)    brcs(f0,t0+2)   brcs(f0+1,t0+2)];

    brcs_weighted = sum([coeff1*sum(term1),coeff2*sum(term2),sum(term3)]);
    nbrcs_scatter = sum(A_eff(nbrcs_doppler_range,nbrcs_delay_range),'all');
    nbrcs_value = brcs_weighted/nbrcs_scatter;

    nbrcs.nbrcs_scatter = nbrcs_scatter;
    nbrcs.nbrcs_value = nbrcs_value;

    % leading edge slope (LES)
    LES_delay_range = t0:t0+num_delay_bins_ddma-1;
    LES_doppler_range = f0-floor(num_doppler_bins_ddma/2):f0+floor(num_doppler_bins_ddma/2);

    LES_region = brcs(LES_doppler_range,LES_delay_range);
    LES_scatter = nbrcs_scatter;
    y_m = sum(LES_region,1);

    tau_chips = [0,0.25,0.5];
    M = length(tau_chips);

    term1 = sum(tau_chips.*y_m);
    term2 = sum(tau_chips)*sum(y_m);
    term3 = sum(tau_chips.*tau_chips);
    term4 = term3*term3;

    alpha_LES = (term1-term2/M)/(term3-term4/M);
    LES_slope = alpha_LES/LES_scatter;

    LES.LES_scatter = LES_scatter;
    LES.LES_slope = LES_slope;

else

    nbrcs.nbrcs_scatter = nan;
    nbrcs.nbrcs_value = nan;

    LES.LES_scatter = nan;
    LES.LES_slope = nan;

end

% trailing edge slope (TES)
brcs_max = max(brcs,[],'all');
[max_doppler_bin,max_delay_bin] = find(brcs==brcs_max,1);

if (max_delay_bin < 37) && (max_doppler_bin < 5) && (max_doppler_bin > 1)

    TES_delay_range = max_delay_bin:max_delay_bin+4;
    TES_doppler_range = max_doppler_bin-1:max_doppler_bin+1;

    TES_region = brcs(TES_doppler_range,TES_delay_range);       % 1 pixel out from peak on trailing edge
    TES_scatter = sum(A_eff(TES_doppler_range,TES_delay_range),'all');

    y_n = sum(TES_region,1);

    tau_chips_TES = [0,0.25,0.5,0.75,1];
    N = length(tau_chips_TES);

    term1_TES = sum(tau_chips_TES.*y_n);
    term2_TES = sum(tau_chips_TES)*sum(y_n);
    term3_TES = sum(tau_chips_TES.*tau_chips_TES);
    term4_TES = term3_TES*term3_TES;

    alpha_TES = (term1_TES-term2_TES/N)/(term3_TES-term4_TES/N);
    TES_slope = alpha_TES/TES_scatter;

    TES.TES_scatter = TES_scatter;
    TES.TES_slope = TES_slope;

else
    TES.TES_scatter = nan;
    TES.TES_slope = nan;

end

    

    
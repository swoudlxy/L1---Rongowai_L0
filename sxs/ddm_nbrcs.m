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

function [nbrcs,LES,TES] = ddm_nbrcs(brcs,A_eff,sp_delay_bin_float,sp_doppler_bin_float)

% define DDMA region size
% TO DO: this region size needs to be refined after getting real data
num_delay_bin_ddma = 5;
num_doppler_bin_ddma = 3;

if (sp_delay_bin_float <= 40) && (sp_delay_bin_float >= 0) && ...
    (sp_doppler_bin_float <= 5) && (sp_doppler_bin_float >= 0)
    
    % sp bin within DDM region
    sp_delay_bin = floor(sp_delay_bin_float);       
    delayOffset_frac = mod(sp_delay_bin_float,sp_delay_bin);

    sp_doppler_bin = floor(sp_doppler_bin_float);
    dopplerOffset_frac = mod(sp_doppler_bin_float,sp_doppler_bin);

    % crop DDMA region
    ddma_delay_range = sp_delay_bin:sp_delay_bin+num_delay_bin_ddma;
    ddma_doppler_range = sp_doppler_bin-ceil(num_doppler_bin_ddma/2):sp_doppler_bin+ceil(num_doppler_bin_ddma/2)+1;

    nbrcs_scatter = A_eff(ddma_doppler_range,ddma_delay_range);     % one more bin more than DDMA dimension

    % crop brcs region - floating SP bin considered
    coeff1 = (1-delayOffset_frac)*(1-dopplerOffset_frac);
    term1 = brcs(sp_doppler_bin-1,sp_delay_bin);

    coeff2 = 1-delayOffset_frac;
    term2 = [   brcs(sp_doppler_bin,sp_delay_bin), ...
                brcs(sp_doppler_bin+1,sp_delay_bin)];

    coeff3 = 1-dopplerOffset_frac;
    term3 = [   brcs(sp_doppler_bin-1,sp_delay_bin+1), ...
                brcs(sp_doppler_bin-1,sp_delay_bin+2), ...
                brcs(sp_doppler_bin-1,sp_delay_bin+3), ...
                brcs(sp_doppler_bin-1,sp_delay_bin+4)];

    coeff4 = delayOffset_frac;
    term4 = [   brcs(sp_doppler_bin,sp_delay_bin+5), ...
                brcs(sp_doppler_bin+1,sp_delay_bin+5)];

    coeff5 = dopplerOffset_frac;
    term5 = [   brcs(sp_doppler_bin+2,sp_delay_bin+1), ...
                brcs(sp_doppler_bin+2,sp_delay_bin+2), ...
                brcs(sp_doppler_bin+2,sp_delay_bin+3), ...
                brcs(sp_doppler_bin+2,sp_delay_bin+4)];

    coeff6 = (1-delayOffset_frac)*dopplerOffset_frac;
    term6 = brcs(sp_doppler_bin+2,sp_delay_bin);

    coeff7 = delayOffset_frac*(1-dopplerOffset_frac);
    term7 = brcs(sp_doppler_bin-1,sp_delay_bin+5);

    coeff8 = delayOffset_frac*dopplerOffset_frac;
    term8 = brcs(sp_doppler_bin-1,sp_delay_bin+5);

    term9 = [   brcs(sp_doppler_bin,sp_delay_bin+1), ...
                brcs(sp_doppler_bin,sp_delay_bin+2), ...
                brcs(sp_doppler_bin,sp_delay_bin+3), ...
                brcs(sp_doppler_bin,sp_delay_bin+4), ...
                brcs(sp_doppler_bin+1,sp_delay_bin+1), ...
                brcs(sp_doppler_bin+1,sp_delay_bin+2), ...
                brcs(sp_doppler_bin+1,sp_delay_bin+3), ...
                brcs(sp_doppler_bin+1,sp_delay_bin+4)];

    brcs_weighted = sum([   coeff1*sum(term1),coeff2*sum(term2),coeff3*sum(term3), ...
                            coeff4*sum(term4),coeff5*sum(term5),coeff6*sum(term6), ...
                            coeff7*sum(term7),coeff8*sum(term8),sum(term9)]);

    nbrcs_value = brcs_weighted/nbrcs_scatter;

    nbrcs.nbrcs_scatter = nbrcs_scatter;
    nbrcs.nbrcs_value = nbrcs_value;

    % leading edge slope (LES)
    LES_delay_range = sp_delay_bin:sp_delay_bin+num_delay_bin_ddma-1;
    LES_doppler_range = sp_doppler_bin-ceil(num_doppler_bin_ddma/2):sp_doppler_bin+ceil(num_doppler_bin_ddma/2);

    LES_region = brcs(LES_doppler_range,LES_delay_range);
    LES_scatter = nbrcs_scatter;
    y_m = sum(LES_region,1);

    tau_chips = [0,0.25,0.5,0.75,1];
    M = length(M);

    term1 = sum(tau_chips.*y_m);
    term2 = sum(tau_chips)*sum(y_m);
    term3 = sum(tau_chips.*tau_chips);
    term4 = term3*term3;

    alpha_LES = (term1-term2/M)/(term3-term4/M);
    LES_slope = alpha_LES/LES_scatter;

    LES.LES_scatter = LES_scatter;
    LES.LES_slope = LES_slope;

    % trailing edge slope (TES)
    TES_delay_range = sp_delay_bin+num_delay_bin_ddma-1:sp_delay_bin+num_delay_bin_ddma;
    
    TES_region = brcs(LES_doppler_range,TES_delay_range);       % 1 pixel out from peak 
    TES_scatter = sum(A_eff(LES_doppler_range,TES_delay_range));

    y_n = sum(TES_region,1);

    tau_chips_TES = [1,1.25];
    N =length(tau_chip_TES);

    term1_TES = sum(tau_chips_TES.*y_n);
    term2_TES = sum(tau_chips_TES)*sum(y_n);
    term3_TES = sum(tau_chips_TES.*tau_chips_TES);
    term4_TES = term3*term3;

    alpha_TES = (term1_TES-term2_TES/N)/(term3_TES-term4_TES/N);
    TES_slope = alpha_TES/TES_scatter;

    TES.TES_scatter = TES_scatter;
    TES.TES_slope = TES_slope;

else

    % sp bin beyond DDM region
    nbrcs.nbrcs_scatter = -9999;
    nbrcs.nbrcs_value = -9999;

    LES.LES_scatter = -9999;
    LES.LES_slope = -9999;

    TES.TES_scatter = -9999;
    TES.TES_slope = -9999;

end
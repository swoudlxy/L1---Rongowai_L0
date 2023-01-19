% this function computes two versions of NBRCS
% floating SP bin location is not considered
% TES and LES are not included in this version

function [nbrcs,nbrcs_scatter] = get_ddm_nbrcs2(brcs,A_eff,sx_bin,flag)

sx_delay_bin = floor(sx_bin(1))+1;
sx_doppler_bin = floor(sx_bin(2))+1;

% case 1: small scattering area
% case 2: large scattering area
switch flag
    case 1
        if sx_delay_bin > 2 && sx_delay_bin <= 40 && ...
                sx_doppler_bin > 1 && sx_doppler_bin < 5

            brcs_ddma = brcs(sx_doppler_bin-1:sx_doppler_bin+1,sx_delay_bin-2:sx_delay_bin);
            A_eff_ddma = A_eff(sx_doppler_bin-1:sx_doppler_bin+1,sx_delay_bin-2:sx_delay_bin);

            brcs_total = sum(brcs_ddma,'all');
            A_eff_total = sum(A_eff_ddma,'all');

            nbrcs = brcs_total/A_eff_total;
            nbrcs_scatter = A_eff_total;

        else
            nbrcs = nan;
            nbrcs_scatter = nan;
        
        end

    case 2
        if sx_delay_bin > 30 && sx_delay_bin <= 40 && ...
                sx_doppler_bin > 1 && sx_doppler_bin < 5

            brcs_ddma = brcs(sx_doppler_bin-1:sx_doppler_bin+1,sx_delay_bin-29:sx_delay_bin);
            A_eff_ddma = A_eff(sx_doppler_bin-1:sx_doppler_bin+1,sx_delay_bin-29:sx_delay_bin);

            brcs_total = sum(brcs_ddma,'all');
            A_eff_total = sum(A_eff_ddma,'all');

            nbrcs = brcs_total/A_eff_total;
            nbrcs_scatter = A_eff_total;

        elseif sx_delay_bin > 1 && sx_delay_bin <= 30 && ...
                sx_doppler_bin > 1 && sx_doppler_bin < 5

            brcs_ddma = brcs(sx_doppler_bin-1:sx_doppler_bin+1,1:sx_delay_bin);
            A_eff_ddma = A_eff(sx_doppler_bin-1:sx_doppler_bin+1,1:sx_delay_bin);

            brcs_total = sum(brcs_ddma,'all');
            A_eff_total = sum(A_eff_ddma,'all');

            nbrcs = brcs_total/A_eff_total;
            nbrcs_scatter = A_eff_total;

        else
            nbrcs = nan;
            nbrcs_scatter = nan;
        
        end



end




    

    
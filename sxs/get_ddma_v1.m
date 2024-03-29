% this function gets the brcs and A_eff within ddma region - SP bin

function [brcs_copol_ddma,brcs_xpol_ddma,A_eff_ddma] = get_ddma_v1(brcs_copol,brcs_xpol,A_eff,sp_delay_row,sp_doppler_col)

delay_intg = floor(sp_delay_row);
delay_frac = sp_delay_row-floor(sp_delay_row);

doppler_intg = round(sp_doppler_col);
doppler_frac = sp_doppler_col-round(sp_doppler_col);

% NBRCS SP bin
if doppler_frac>=0

    if doppler_intg<=4

        brcs_copol_ddma =  (1-doppler_frac)*(1-delay_frac)*brcs_copol(doppler_intg,delay_intg)+ ...
                           (1-doppler_frac)*delay_frac*brcs_copol(doppler_intg,delay_intg+1)+ ...
                           doppler_frac*(1-delay_frac)*brcs_copol(doppler_intg+1,delay_intg)+ ...
                           doppler_frac*delay_frac*brcs_copol(doppler_intg+1,delay_intg+1);

        brcs_xpol_ddma =   (1-doppler_frac)*(1-delay_frac)*brcs_xpol(doppler_intg,delay_intg)+ ...
                           (1-doppler_frac)*delay_frac*brcs_xpol(doppler_intg,delay_intg+1)+ ...
                           doppler_frac*(1-delay_frac)*brcs_xpol(doppler_intg+1,delay_intg)+ ...
                           doppler_frac*delay_frac*brcs_xpol(doppler_intg+1,delay_intg+1);

        A_eff_ddma      =  (1-doppler_frac)*(1-delay_frac)*A_eff(doppler_intg,delay_intg)+ ...
                           (1-doppler_frac)*delay_frac*A_eff(doppler_intg,delay_intg+1)+ ...
                           doppler_frac*(1-delay_frac)*A_eff(doppler_intg+1,delay_intg)+ ...
                           doppler_frac*delay_frac*A_eff(doppler_intg+1,delay_intg+1);

    elseif doppler_intg>4

        brcs_copol_ddma =  (1-delay_frac)*brcs_copol(doppler_intg,delay_intg)+ ...
                           delay_frac*brcs_copol(doppler_intg,delay_intg+1);

        brcs_xpol_ddma =   (1-delay_frac)*brcs_xpol(doppler_intg,delay_intg)+ ...
                           delay_frac*brcs_xpol(doppler_intg,delay_intg+1);

        A_eff_ddma      =  (1-delay_frac)*A_eff(doppler_intg,delay_intg)+ ...
                           delay_frac*A_eff(doppler_intg,delay_intg+1);

    end

elseif doppler_frac<0

    if doppler_intg>=2

        brcs_copol_ddma =  (1-abs(doppler_frac))*(1-delay_frac)*brcs_copol(doppler_intg,delay_intg)+ ...
                           (1-abs(doppler_frac))*delay_frac*brcs_copol(doppler_intg,delay_intg+1)+ ...
                           abs(doppler_frac)*(1-delay_frac)*brcs_copol(doppler_intg-1,delay_intg)+ ...
                           abs(doppler_frac)*delay_frac*brcs_copol(doppler_intg-1,delay_intg+1);

        brcs_xpol_ddma =   (1-abs(doppler_frac))*(1-delay_frac)*brcs_xpol(doppler_intg,delay_intg)+ ...
                           (1-abs(doppler_frac))*delay_frac*brcs_xpol(doppler_intg,delay_intg+1)+ ...
                           abs(doppler_frac)*(1-delay_frac)*brcs_xpol(doppler_intg-1,delay_intg)+ ...
                           abs(doppler_frac)*delay_frac*brcs_xpol(doppler_intg-1,delay_intg+1);

        A_eff_ddma      =  (1-abs(doppler_frac))*(1-delay_frac)*A_eff(doppler_intg,delay_intg)+ ...
                           (1-abs(doppler_frac))*delay_frac*A_eff(doppler_intg,delay_intg+1)+ ...
                           abs(doppler_frac)*(1-delay_frac)*A_eff(doppler_intg-1,delay_intg)+ ...
                           abs(doppler_frac)*delay_frac*A_eff(doppler_intg-1,delay_intg+1);

    elseif doppler_intg<2

        brcs_copol_ddma = (1-delay_frac)*brcs_copol(doppler_intg,delay_intg)+ ...
                           delay_frac*brcs_copol(doppler_intg,delay_intg+1);

        brcs_xpol_ddma =  (1-delay_frac)*brcs_xpol(doppler_intg,delay_intg)+ ...
                           delay_frac*brcs_xpol(doppler_intg,delay_intg+1);

        A_eff_ddma      = (1-delay_frac)*A_eff(doppler_intg,delay_intg)+ ...
                           delay_frac*A_eff(doppler_intg,delay_intg+1);

    end

end
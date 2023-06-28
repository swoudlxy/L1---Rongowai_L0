function ddm = rongowaiWAF(doppler_center,delay_center,dim,delay_axis_chips)
% dim = [DelayBins DopplerBins]
cohT = 1e-3;
tauChip_s = 1/1.023e6;
delaybins = dim(2)-1;
dopplerbins = dim(1)-1;
centerDoppler = doppler_center-1;
centerDelay = delay_center-1;
resChips = 0.2552/2;
resDoppler = 250;%500;           % Rongowai Doppler Res should be 250

for l=0:delaybins
    dtau_s = delay_axis_chips(l+1)*tauChip_s;
    for k=0:dopplerbins
        dfreq_Hz = (k-centerDoppler)*resDoppler;
        if abs(dtau_s)<=(tauChip_s*(1+tauChip_s/cohT))
            lambda = 1-abs(dtau_s)/tauChip_s;
        else
            lambda = -tauChip_s/cohT;
        end
        x = dfreq_Hz*pi*cohT;
        rads = -pi*dfreq_Hz*cohT;
        if x==0
            s = 1;
        else
            s = (sin(x)/x)*exp(1i*rads);
        end
        ddm(l+1,k+1) = lambda*s;
    end
end
ddm = abs(ddm).'; 
ddm = ddm.*ddm;

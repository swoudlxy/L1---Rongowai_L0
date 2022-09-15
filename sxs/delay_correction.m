% this function correct the input code phase to a value between 0 and
% a defined value P, P = 1023 for CPS L1 C/A

function delay_chips_out = delay_correction(delay_chips_in,P)

temp = delay_chips_in;

if temp < 0
    while temp < 0
        temp = temp+P;
    end  

elseif temp > 1023
    while temp >1023
        temp = temp-P;
    end

end

delay_chips_out = temp;
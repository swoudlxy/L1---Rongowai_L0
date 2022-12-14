% this function derives the confidence flag for the computed SP and also
% the zenith code phase directly tracked by the NGRx

function [specular_bin,zenith_code_phase,confidence_flag] = get_specular_bin(tx,rx,sx,ddm)

c = 299792458;

tx_pos_xyz = tx.tx_pos_xyz;
tx_vel_xyz = tx.tx_vel_xyz;

rx_pos_xyz = rx.rx_pos_xyz;
rx_vel_xyz = rx.rx_vel_xyz;
rx_clk_drift = rx.rx_clk_drift;

sx_pos_xyz = sx.sx_pos_xyz;
sx_d_snell = sx.sx_d_snell;
dist_to_coast = sx.dist_to_coast;

raw_counts = ddm.raw_counts;
delay_resolution = ddm.delay_resolution;
delay_center_chips = ddm.delay_center_chips;
delay_center_bin = ddm.delay_center_bin;

doppler_resolution = ddm.doppler_resolution;
doppler_center_hz = ddm.doppler_center_hz;
doppler_center_bin = ddm.doppler_center_bin;

add_range_to_sp = ddm.add_range_to_sp;
snr_db = ddm.snr_db;

% derive zenith code phase
add_range_to_sp_chips = meter2chips(add_range_to_sp);
zenith_code_phase1 = delay_center_chips+add_range_to_sp_chips;
zenith_code_phase = delay_correction(zenith_code_phase1,1023);

% derive precise SP bin location
[pixel_delay_chips,pixel_doppler_hz,pixel_add_range_to_sp_chips] = deldop(tx_pos_xyz,rx_pos_xyz, ...
    tx_vel_xyz,rx_vel_xyz,sx_pos_xyz);

delay_error = add_range_to_sp_chips-pixel_add_range_to_sp_chips;
sp_delay_row = delay_center_bin+delay_error/delay_resolution;           % cygnss is using "-"

doppler_clk = rx_clk_drift/c;
pixel_doppler_hz = pixel_doppler_hz+doppler_clk;

doppler_error = doppler_center_hz-pixel_doppler_hz;
sp_dopp_col = doppler_center_bin-doppler_error/doppler_resolution;

sp_delay_error = delay_center_chips-pixel_delay_chips;
sp_dopp_error = doppler_center_hz-pixel_doppler_hz;

% derive confidence flag
if dist_to_coast < 0
    confidence_flag = 3;

else
    peak_counts = max(raw_counts,[],'all');
    [doppler_max_bin,delay_max_bin] = find(raw_counts == peak_counts,1);

    delay_max = delay_center_chips+(delay_max_bin-delay_center_bin-1)*delay_resolution;
    delay_sp = zenith_code_phase1-pixel_add_range_to_sp_chips;
    delay_diff = abs(delay_sp-delay_max);

    doppler_max = doppler_center_hz+(doppler_max_bin-doppler_center_bin-1)*doppler_resolution;
    doppler_diff = abs(pixel_doppler_hz-doppler_max);

    delay_doppler_snell = (delay_diff<2.5) && (doppler_diff<200) && (sx_d_snell<2);

    if snr_db>=2 && ~delay_doppler_snell
        confidence_flag = 0;
    elseif snr_db<2 && ~delay_doppler_snell
        confidence_flag = 1;
    elseif snr_db<2 && delay_doppler_snell
        confidence_flag = 3;
    elseif snr_db>=2 && delay_doppler_snell
        confidence_flag = 4;
    else
        confidence_flag = nan;
    end

end

specular_bin = [sp_delay_row,sp_dopp_col,sp_delay_error,sp_dopp_error];


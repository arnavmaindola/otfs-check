%% main_csi.m - BER vs velocity for several CSI qualities (OTFS, OFDM one-tap, OFDM full MMSE)
% Receiver uses H_est = H + E, E complex Gaussian i.i.d. with per-entry variance
% sigma_e^2 = 10^(dB/10) * mean(|H(:)|^2) (same for the matrix each detector inverts;
% the one-tap detector uses diag(H_est)). The true H and the noise drive the
% transmission. Same seeds/draw order as sim_link, so the 'perfect' curve equals
% the perfect-CSI A2 result (E is drawn AFTER the channel/bits/noise draws).
clear; clc; close all;

quick = true;                          % true: 20 frames, 3 velocities (smoke test)
NUM_FRAMES = 300;                      % Monte-Carlo frames per point
N_GRID = 8;                            % Doppler bins N
csi_db = [-Inf -30 -20 -15 -10];       % CSI error power (dB): OTFS/full-MMSE rel. to mean |H|^2 entry; one-tap rel. to subcarrier gain power; -Inf = perfect

num_frames = NUM_FRAMES; vel_list = 0:100:700; snr_db = 20; seed = 1000;
if quick, num_frames = 20; vel_list = [0 300 700]; end

params.M = 64;  params.N = N_GRID;  params.delta_f = 15e3;  params.fc = 4e9;
params.velocity_kmh = 120; params.snr_db = snr_db; params.fading_model = 'rayleigh';
M = params.M; N = params.N; MN = M*N; nbits = 2*MN;
nl = numel(csi_db); nvel = numel(vel_list);
err = zeros(nvel, nl, 3);              % schemes: 1 OTFS, 2 OFDM one-tap, 3 OFDM full MMSE

tic;
for iv = 1:nvel
    pA = params; pA.velocity_kmh = vel_list(iv);
    for f = 1:num_frames
        rng(seed + f);
        chan = channel_state(pA);
        p = pA; p.cp_len = chan.max_delay;
        bits = randi([0 1], nbits, 1);
        tx_o = otfs_transmitter(p, bits);
        tx_f = ofdm_transmitter(p, bits);
        [y_o, H_o, nv] = otfs_channel(tx_o, p, chan);
        [y_f, H_f, ~]  = ofdm_channel(tx_f, p, chan);
        % error draws after all perfect-path randomness; shared across levels (nested)
        E_o = (randn(MN) + 1i*randn(MN)) / sqrt(2);
        E_f = (randn(MN) + 1i*randn(MN)) / sqrt(2);
        E_d = (randn(MN,1) + 1i*randn(MN,1)) / sqrt(2);   % one-tap per-subcarrier error
        pw_o = mean(abs(H_o(:)).^2);  pw_f = mean(abs(H_f(:)).^2);
        pw_d = mean(abs(diag(H_f)).^2);                   % per-subcarrier gain power
        for l = 1:nl
            if isinf(csi_db(l)), s = 0; else, s = 10^(csi_db(l)/20); end
            He_o = H_o + s*sqrt(pw_o)*E_o;
            He_f = H_f + s*sqrt(pw_f)*E_f;
            He_1 = H_f + diag(s*sqrt(pw_d)*E_d);
            [~, ~, ~, e] = otfs_receiver(y_o, He_o, nv, bits, p);
            err(iv, l, 1) = err(iv, l, 1) + e;
            [~, ~, ~, e] = ofdm_receiver(y_f, He_1, nv, bits, p, 'onetap');
            err(iv, l, 2) = err(iv, l, 2) + e;
            [~, ~, ~, e] = ofdm_receiver(y_f, He_f, nv, bits, p, 'mmse');
            err(iv, l, 3) = err(iv, l, 3) + e;
        end
    end
end
ber = err / (nbits * num_frames);

tag = sprintf('CSI_N%d%s', N, ternary(quick, '_quick', ''));
scheme = {'OTFS (MMSE)', 'OFDM (one-tap)', 'OFDM (full MMSE)'};
lbl = arrayfun(@(d) ternary(isinf(d), 'perfect CSI', sprintf('%d dB error', d)), csi_db, 'UniformOutput', false);
capt = sprintf('M=%d, N=%d, f_c=%g GHz, SNR=%d dB, %d frames', M, N, params.fc/1e9, snr_db, num_frames);
cols = lines(nl); cols(1,:) = [0 0 0];
st = arrayfun(@(l) {cols(l,:), '-o'}, 1:nl, 'UniformOutput', false);

diary([tag '_log.txt']); diary on;
fprintf('%s; hollow = 0 errors\n', capt);
for k = 1:3
    fprintf('\n%s: BER vs velocity\n km/h  %s\n', scheme{k}, sprintf('%-14s', lbl{:}));
    for iv = 1:nvel
        fprintf('%4d  %s\n', vel_list(iv), sprintf('%-14.2e', ber(iv,:,k)));
    end
    fig = figure;
    h = plot_ber(vel_list(:), ber(:,:,k), nbits*num_frames, 1e-5, st);
    xlabel('Velocity (km/h)'); ylabel('BER');
    title({[scheme{k} ': BER vs velocity, CSI quality'], capt}, 'FontSize', 8);
    legend(h, lbl, 'Location', 'northwest');
    saveas(fig, sprintf('%s_%d.png', tag, k));
end
save([tag '_results.mat'], 'ber', 'vel_list', 'csi_db', 'num_frames', 'M', 'N', 'snr_db');
fprintf('\nTotal time: %.1f s\n', toc);
diary off;

function s = ternary(c, a, b)
    if c, s = a; else, s = b; end
end

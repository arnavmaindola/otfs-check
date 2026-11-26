%% main_csi_final.m - EXTRA (beyond synopsis): imperfect CSI at M=32, N=16
% Receiver CSI: H_est = H + s*sqrt(mean|H|^2)*E, E ~ CN(0,1) i.i.d. (one-tap: error on
% the diagonal, relative to subcarrier-gain power). True H and noise drive the link.
%   mode = 'A' : fixed error levels csi_db, BER vs velocity 0:50:500 @ 20 dB
%   mode = 'B' : SNR-tied error, NMSE = 1/(Lp*SNR_lin) (pilot LS-style); BER vs SNR
%                at 120 and 500 km/h, Lp in Lp_list, perfect CSI as reference
% Same bits/noise draws across CSI levels. Perfect-CSI curve uses this script's draw
% order (matches sim_link only approximately for mode B, exactly in distribution).
% Runtime (1000 frames): A ~ 60-90 min, B ~ 60-100 min; run quick = true first (~2 min).
clear; clc; close all;

mode = 'A';                            % 'A' or 'B'
quick = true;
NUM_FRAMES = 1000;
N_GRID = 16;  M_GRID = 32;
csi_db  = [-Inf -30 -20 -15 -10];      % mode A levels (dB); -Inf = perfect
Lp_list = [1 10];                      % mode B pilot-gain factors (NMSE = 1/(Lp*SNR))
snr_A = 20; vel_A = 0:50:500; snr_B = 0:5:30; vel_B = [120 500];

num_frames = NUM_FRAMES; if quick, num_frames = 20; vel_A = [0 250 500]; snr_B = 0:10:30; end
params.M = M_GRID; params.N = N_GRID; params.delta_f = 15e3; params.fc = 4e9;
params.fading_model = 'rayleigh';
M = M_GRID; N = N_GRID; MN = M*N; nb = 2*MN; seed = 1000;

if mode == 'A'
    sweeps = vel_A(:); nl = numel(csi_db); lbl = arrayfun(@(d) ternary(isinf(d),'perfect CSI',sprintf('%d dB error',d)), csi_db, 'UniformOutput', false);
    xl = 'Velocity (km/h)'; fixed = sprintf('SNR=%d dB', snr_A);
else
    nl = 1 + numel(Lp_list); lbl = [{'perfect CSI'}, arrayfun(@(l) sprintf('NMSE=1/(%d\\cdotSNR)', l), Lp_list, 'UniformOutput', false)];
    xl = 'SNR (dB)';
end
scheme = {'OTFS (MMSE)', 'OFDM (one-tap)', 'OFDM (full MMSE)'};
tag = sprintf('CSI%s_M%d_N%d%s', mode, M, N, ternary(quick, '_quick', ''));
diary([tag '_log.txt']); diary on; tic;

if mode == 'A', groups = 1; else, groups = numel(vel_B); end
res = cell(groups,1);
for g = 1:groups
    if mode == 'A', xs = vel_A(:); else, xs = snr_B(:); end
    err = zeros(numel(xs), nl, 3);
    for ix = 1:numel(xs)
        p = params;
        if mode == 'A', p.velocity_kmh = xs(ix); p.snr_db = snr_A; snr_now = snr_A;
        else, p.velocity_kmh = vel_B(g); p.snr_db = xs(ix); snr_now = xs(ix); end
        for f = 1:num_frames
            rng(seed + f);
            chan = channel_state(p); p.cp_len = chan.max_delay;
            bits = randi([0 1], nb, 1);
            [y_o, H_o, nv] = otfs_channel(otfs_transmitter(p, bits), p, chan);
            [y_f, H_f, ~]  = ofdm_channel(ofdm_transmitter(p, bits), p, chan);
            E_o = (randn(MN)+1i*randn(MN))/sqrt(2); E_f = (randn(MN)+1i*randn(MN))/sqrt(2);
            E_d = (randn(MN,1)+1i*randn(MN,1))/sqrt(2);
            for l = 1:nl
                if mode == 'A'
                    if isinf(csi_db(l)), s = 0; else, s = 10^(csi_db(l)/20); end
                else
                    if l == 1, s = 0; else, s = sqrt(1/(Lp_list(l-1)*10^(snr_now/10))); end
                end
                err(ix,l,:) = squeeze(err(ix,l,:)) + csi_eval(y_o,H_o,y_f,H_f,nv,bits,p,s,s,s,E_o,E_f,E_d)';
            end
        end
    end
    res{g} = err / (nb*num_frames);
end

nbits_total = nb*num_frames; cols = lines(nl); cols(1,:) = [0 0 0];
st = arrayfun(@(l) {cols(l,:), '-o'}, 1:nl, 'UniformOutput', false);
for g = 1:groups
    if mode == 'A', xs = vel_A(:); fx = sprintf('SNR %d dB', snr_A); else, xs = snr_B(:); fx = sprintf('%d km/h', vel_B(g)); end
    capt = sprintf('M=%d N=%d fc=4GHz %s %d frames', M, N, fx, num_frames);
    for k = 1:3
        b = res{g}(:,:,k);
        fprintf('\n%s, %s\n x     %s\n', scheme{k}, capt, sprintf('%-16s', lbl{:}));
        for i = 1:numel(xs), fprintf('%4d  %s\n', xs(i), sprintf('%-16.2e', b(i,:))); end
        fig = figure; h = plot_ber(xs, b, nbits_total, 5e-7, st);
        xlabel(xl); ylabel('BER'); legend(h, lbl, 'Location', 'best', 'FontSize', 8);
        title({[scheme{k} ' - imperfect CSI'], capt}, 'FontSize', 8);
        saveas(fig, sprintf('%s_g%d_s%d.png', tag, g, k));
    end
end
save([tag '_results.mat'], 'res', 'mode', 'num_frames', 'csi_db', 'Lp_list', 'vel_A', 'snr_B', 'vel_B');
fprintf('\nTotal time: %.1f s\n', toc); diary off;

function s = ternary(c, a, b)
    if c, s = a; else, s = b; end
end

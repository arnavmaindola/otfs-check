%% main_partB_final.m - FINAL Part B: OTFS under 5 fading models (M=32, N=16)
% B1: BER vs SNR @ 120 km/h     B2: BER vs velocity 0:50:500 km/h @ 15 dB
% Models and parameters (all E|h_p|^2 = EVA power p_p):
%   Rayleigh
%   Rician      K = 6 dB (3.98), LOS ray on tap 1 only
%   Nakagami-m  m = 2, every tap
%   alpha-mu    alpha = 1.2, mu = 0.8 (severe), every tap
%   kappa-mu    kappa = K of the Rician = 3.98 (6 dB), mu = 2, every tap
%     (the channel_state default kappa=3 is overridden so kappa equals the Rician K)
% CRN pairing: every model is built from the SAME per-frame variates, bits and noise;
% error bars (95% CI) and paired differences vs Rayleigh quantify resolvability.
% Runtime at 1000 frames: ~35-50 min (5 models x 18 points); run quick = true first (~1 min).
clear; clc; close all;

quick = true;                          % true: 40 frames, coarse grids
NUM_FRAMES = 1000;
N_GRID = 16;  M_GRID = 32;
num_frames = NUM_FRAMES; if quick, num_frames = 40; end
snr_list = 0:5:30; vel_list = 0:50:500; v_snr = 120; snr_for_vel = 15;
if quick, snr_list = [0 15 30]; vel_list = [0 250 500]; end
seed = 5000;

K = 10^(6/10);
base.M = M_GRID; base.N = N_GRID; base.delta_f = 15e3; base.fc = 4e9; base.snr_db = 10;
cfg = {'Rayleigh',                       struct('fading_model','rayleigh');
       'Rician (K=6 dB, tap 1)',         struct('fading_model','rician');
       'Nakagami-m (m=2)',               struct('fading_model','nakagami-m');
       '\alpha-\mu (\alpha=1.2, \mu=0.8)', struct('fading_model','alpha-mu','alpha',1.2,'mu',0.8);
       sprintf('\\kappa-\\mu (\\kappa=%.2f=K_{Rice}, \\mu=2)', K), struct('fading_model','kappa-mu','kappa',K,'mu',2)};
nc = size(cfg,1); P = 9;
plain = regexprep(cfg(:,1), '[\\{}]', '');

pts = [snr_list(:), v_snr*ones(numel(snr_list),1);
       snr_for_vel*ones(numel(vel_list),1), vel_list(:)];
npt = size(pts,1); nsn = numel(snr_list);
nb = 2*base.M*base.N;
fb = zeros(num_frames, npt, nc);

tag0 = sprintf('B_M%d_N%d%s', M_GRID, N_GRID, ternary(quick, '_quick', ''));
diary([tag0 '_log.txt']); diary on; tic;
capt = sprintf('OTFS MMSE, CRN-paired, M=%d, N=%d, f_c=%g GHz, %d frames', base.M, base.N, base.fc/1e9, num_frames);
fprintf('%s\n', capt);

for f = 1:num_frames
    rng(seed + f);
    V.ang = rand(P,1); V.ph = rand(P,1); V.u = min(max(rand(P,1),1e-12),1-1e-12); V.n = randn(8,P);
    bits = randi([0 1], nb, 1);
    tx = otfs_transmitter(base, bits);
    for c = 1:nc
        p = base; fn = fieldnames(cfg{c,2});
        for q = 1:numel(fn), p.(fn{q}) = cfg{c,2}.(fn{q}); end
        p.crn = V;
        for k = 1:npt
            p.snr_db = pts(k,1); p.velocity_kmh = pts(k,2);
            chan = channel_state(p); p.cp_len = chan.max_delay;
            rng(seed + 100000*f + k);                       % identical noise across models
            [y, H, nv] = otfs_channel(tx, p, chan);
            [~,~,~,e] = otfs_receiver(y, H, nv, bits, p);
            fb(f,k,c) = e / nb;
        end
    end
end

ylo = 1e-5; nbits_total = nb*num_frames;
mean_b = squeeze(mean(fb,1));
se_b = squeeze(std(fb,0,1)) / sqrt(num_frames);
for part = 1:2
    if part == 1, rows = 1:nsn; xs = snr_list(:); xl = 'SNR (dB)'; nm = 'B1'; ttl = sprintf('BER vs SNR (%d km/h)', v_snr); loc = 'southwest';
    else, rows = nsn+(1:numel(vel_list)); xs = vel_list(:); xl = 'Velocity (km/h)'; nm = 'B2'; ttl = sprintf('BER vs velocity (SNR %d dB)', snr_for_vel); loc = 'best'; end
    fprintf('\n[%s] OTFS %s\n%5s %s\n', nm, ttl, 'x', sprintf('| %-24s', plain{:}));
    for r = 1:numel(rows)
        cells = arrayfun(@(c) sprintf('%.3e +-%.1e', mean_b(rows(r),c), se_b(rows(r),c)), 1:nc, 'UniformOutput', false);
        fprintf('%5d %s\n', xs(r), sprintf('| %-24s', cells{:}));
    end
    fprintf('Paired difference vs Rayleigh (mean diff, z):\n');
    for a = 2:nc
        d = squeeze(fb(:,rows,a) - fb(:,rows,1));
        md = mean(d,1); sd = std(d,0,1)/sqrt(num_frames); z = md./max(sd,eps);
        fprintf('  %-34s %s\n', plain{a}, sprintf('%d:%+.1e(z%+.1f) ', [xs.'; md; z]));
    end
    fig = figure; cols = lines(nc);
    st = arrayfun(@(k) {cols(k,:), '-o'}, 1:nc, 'UniformOutput', false);
    h = plot_ber(xs, mean_b(rows,:), nbits_total, ylo, st);
    for k = 1:nc
        errorbar(xs, mean_b(rows,k), 2*se_b(rows,k), 'Color', cols(k,:), 'LineStyle', 'none', 'HandleVisibility', 'off');
    end
    ylim([ylo 1]); xlabel(xl); ylabel('BER'); legend(h, cfg(:,1), 'Location', loc, 'FontSize', 8);
    title({sprintf('%s: OTFS %s (bars = 95%% CI)', nm, ttl), capt}, 'FontSize', 8);
    saveas(fig, sprintf('%s_%s.png', tag0, nm));
end
save([tag0 '_results.mat'], 'fb', 'pts', 'cfg', 'num_frames', 'snr_list', 'vel_list');
fprintf('\nTotal time: %.1f s\n', toc);
diary off;

function s = ternary(c, a, b)
    if c, s = a; else, s = b; end
end

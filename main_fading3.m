%% main_fading3.m - OTFS fading-model study with COMMON RANDOM NUMBERS (paired comparison)
% Per frame, ONE set of variates (Doppler angles, tap phases, 8xP normals, uniforms, bits,
% noise seed) is shared by every fading model; each model builds its gains from the same
% variates (channel_state 'crn' hook; default behaviour of channel_state is unchanged).
% Differences between models are then paired per frame, so their standard errors are far
% smaller than those of the independent BERs. kappa-mu mu=1 and 'rician-all-taps' have the
% same distribution (not identical samples: kappa-mu uses a fixed LOS phase + independent
% uniform phase, Rician a random LOS phase), which is a validity check (expect |z|<2).
% Theory: kappa-mu AF = (1+2k)/(mu(1+k)^2), k=3.98: mu=1 0.36, 2 0.18, 4 0.09 -> BER falls
% with mu; all alpha-mu sets below are worse than Rayleigh. Judge at SNR >= 15 dB.
% OTFS + MMSE, perfect CSI. Runtime ~ 0.4 s per frame (10 models x 15 points).
clear; clc; close all;

quick = true;                          % true: 40 frames, coarse grids
NUM_FRAMES = 1000;                     % Monte-Carlo frames (CRN-paired)
N_GRID = 8;
num_frames = NUM_FRAMES; if quick, num_frames = 40; end
snr_list = 0:5:30; vel_list = 0:100:700; v_snr = 120; snr_for_vel = 15;
if quick, snr_list = [0 15 30]; vel_list = [0 300 700]; end
seed = 5000;

K = 10^(6/10);
base.M = 64; base.N = N_GRID; base.delta_f = 15e3; base.fc = 4e9; base.snr_db = 10;
cfg = {'Rayleigh',                   struct('fading_model','rayleigh');
       'Rician K=6dB (tap 1 only)',  struct('fading_model','rician');
       'Rician K=6dB (all taps)',    struct('fading_model','rician-all-taps');
       'kappa-mu k=3.98, mu=1',      struct('fading_model','kappa-mu','kappa',K,'mu',1);
       'kappa-mu k=3.98, mu=2',      struct('fading_model','kappa-mu','kappa',K,'mu',2);
       'kappa-mu k=3.98, mu=4',      struct('fading_model','kappa-mu','kappa',K,'mu',4);
       'alpha-mu (1,0.5)',           struct('fading_model','alpha-mu','alpha',1,  'mu',0.5);
       'alpha-mu (2,0.5)',           struct('fading_model','alpha-mu','alpha',2,  'mu',0.5);
       'alpha-mu (1.2,0.8)',         struct('fading_model','alpha-mu','alpha',1.2,'mu',0.8);
       'alpha-mu (1,0.3)',           struct('fading_model','alpha-mu','alpha',1,  'mu',0.3)};
nc = size(cfg,1); P = 9;
studies = {'a_kappa_mu', [1 2 3 4 5 6]; 'b_alpha_mu', [1 7 8 9 10]};
% paired contrasts (index pairs a-b) reported per study
contr = {'a_kappa_mu', [2 1; 3 1; 4 1; 4 3; 5 1; 6 1; 5 4; 6 5];
         'b_alpha_mu', [7 1; 8 1; 9 1; 10 1]};

pts = [snr_list(:), v_snr*ones(numel(snr_list),1);          % SNR sweep @ v_snr
       snr_for_vel*ones(numel(vel_list),1), vel_list(:)];    % velocity sweep @ snr_for_vel
npt = size(pts,1); nsn = numel(snr_list);
nb = 2*base.M*base.N; MN = base.M*base.N;
fb = zeros(num_frames, npt, nc);

tag0 = sprintf('FAD3_N%d%s', N_GRID, ternary(quick, '_quick', ''));
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
mean_b = squeeze(mean(fb,1));                       % npt x nc
se_b = squeeze(std(fb,0,1)) / sqrt(num_frames);     % unpaired-style SE of each curve
R = struct();
for s = 1:size(studies,1)
    nm = studies{s,1}; idx = studies{s,2};
    lbl = cfg(idx,1); plain = lbl;
    ct = contr{s,2};
    for part = 1:2
        if part == 1, rows = 1:nsn; xs = snr_list(:); xl = 'SNR (dB)'; ttl = sprintf('BER vs SNR (%d km/h)', v_snr); loc = 'southwest';
        else, rows = nsn+(1:numel(vel_list)); xs = vel_list(:); xl = 'Velocity (km/h)'; ttl = sprintf('BER vs velocity (SNR %d dB)', snr_for_vel); loc = 'best'; end
        fprintf('\n[%s] %s\n%5s %s\n', nm, ttl, 'x', sprintf('| %-24s', plain{:}));
        for r = 1:numel(rows)
            cells = arrayfun(@(c) sprintf('%.3e +-%.1e', mean_b(rows(r),c), se_b(rows(r),c)), idx, 'UniformOutput', false);
            fprintf('%5d %s\n', xs(r), sprintf('| %-24s', cells{:}));
        end
        fprintf('Paired differences (mean diff +- paired SE, z):\n');
        for j = 1:size(ct,1)
            a = ct(j,1); b = ct(j,2);
            d = squeeze(fb(:,rows,a) - fb(:,rows,b));      % frames x points
            md = mean(d,1); sd = std(d,0,1)/sqrt(num_frames); z = md./max(sd,eps);
            fprintf('  %-26s - %-24s: %s\n', cfg{a,1}, cfg{b,1}, ...
                sprintf('%d:%+.1e(z%+.1f) ', [xs.'; md; z]));
        end
        fig = figure; cols = lines(numel(idx));
        st = arrayfun(@(k) {cols(k,:), '-o'}, 1:numel(idx), 'UniformOutput', false);
        h = plot_ber(xs, mean_b(rows, idx), nbits_total, ylo, st);
        for k = 1:numel(idx)
            errorbar(xs, mean_b(rows,idx(k)), 2*se_b(rows,idx(k)), 'Color', cols(k,:), 'LineStyle', 'none', 'HandleVisibility', 'off');
        end
        ylim([ylo 1]); xlabel(xl); ylabel('BER'); legend(h, lbl, 'Location', loc, 'FontSize', 8);
        title({sprintf('Fading study %s: %s (bars = 95%% CI)', nm, ttl), capt}, 'FontSize', 8);
        saveas(fig, sprintf('%s_%s_%d.png', tag0, nm, part));
    end
end
save([tag0 '_results.mat'], 'fb', 'pts', 'cfg', 'num_frames', 'snr_list', 'vel_list');
fprintf('\nTotal time: %.1f s\n', toc);
diary off;

function s = ternary(c, a, b)
    if c, s = a; else, s = b; end
end

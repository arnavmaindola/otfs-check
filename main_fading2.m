%% main_fading2.m - OTFS fading-model study with deliberately chosen parameters
% (a) kappa-mu with kappa = K of Rician (K = 6 dB) and mu in {1,2,4}, vs Rician and Rayleigh
% (b) alpha-mu with poor-fading parameter sets, vs Rayleigh
% BER vs SNR @ v_snr km/h and BER vs velocity @ snr_for_vel dB, OTFS + MMSE, perfect CSI.
%
% Definitions: all models have E|h_p|^2 = EVA power p_p.
%  'rician'          : LOS ray on tap 1 ONLY (K = 6 dB of total power)  [existing model]
%  'rician-all-taps' : K = 6 dB on EVERY tap (== kappa-mu with mu = 1, like-for-like)
%  'kappa-mu'        : on EVERY tap; kappa = K, mu = number of clusters
%  'alpha-mu'        : every tap, power-normalised, random phase
% Theory (amount of fading AF): kappa-mu AF = (1+2k)/(mu(1+k)^2) -> falls as mu grows
%  (k=3.98: mu=1 0.36, mu=2 0.18, mu=4 0.09), so expected BER order (high SNR):
%  Rayleigh > rician-all-taps = kappa-mu(mu=1) > mu=2 > mu=4. 'rician' (tap-1 only) lies
%  between Rayleigh and rician-all-taps. alpha-mu: AF grows as mu falls / alpha changes,
%  all sets should be worse than Rayleigh at high SNR (smaller mu = worse).
% At very low SNR (0-5 dB) BER curves can cross (BER concave in SINR); judge at >= 15 dB.
clear; clc; close all;

quick = true;                          % true: 20 frames, coarse grids (smoke test)
NUM_FRAMES = 300;                      % Monte-Carlo frames per point
N_GRID = 8;
num_frames = NUM_FRAMES; if quick, num_frames = 20; end

K_lin = 10^(6/10);                     % kappa = K of the Rician (6 dB)
snr_list = 0:5:30; vel_list = 0:100:700; v_snr = 120; snr_for_vel = 15;
if quick, snr_list = [0 15 30]; vel_list = [0 300 700]; end

params.M = 64;  params.N = N_GRID;  params.delta_f = 15e3;  params.fc = 4e9;
params.velocity_kmh = v_snr; params.snr_db = 10; params.fading_model = 'rayleigh';

% each study: list of {label, fields to set}
cfg_ref = {'Rayleigh',                  struct('fading_model','rayleigh');
           'Rician K=6dB (tap 1 only)', struct('fading_model','rician');
           'Rician K=6dB (all taps)',   struct('fading_model','rician-all-taps')};
cfg_a = [cfg_ref;
    {sprintf('\\kappa-\\mu \\kappa=%.2f, \\mu=1', K_lin), struct('fading_model','kappa-mu','kappa',K_lin,'mu',1);
     sprintf('\\kappa-\\mu \\kappa=%.2f, \\mu=2', K_lin), struct('fading_model','kappa-mu','kappa',K_lin,'mu',2);
     sprintf('\\kappa-\\mu \\kappa=%.2f, \\mu=4', K_lin), struct('fading_model','kappa-mu','kappa',K_lin,'mu',4)}];
cfg_b = {'Rayleigh', struct('fading_model','rayleigh');
         '\alpha-\mu (1,0.5)',     struct('fading_model','alpha-mu','alpha',1,  'mu',0.5);
         '\alpha-\mu (2,0.5)',     struct('fading_model','alpha-mu','alpha',2,  'mu',0.5);
         '\alpha-\mu (1.2,0.8)',   struct('fading_model','alpha-mu','alpha',1.2,'mu',0.8);
         '\alpha-\mu (1,0.3)',     struct('fading_model','alpha-mu','alpha',1,  'mu',0.3)};
studies = {'a_kappa_mu', cfg_a; 'b_alpha_mu', cfg_b};

nbits_total = 2*params.M*params.N*num_frames;
capt = sprintf('OTFS, MMSE, M=%d, N=%d, f_c=%g GHz, %d frames', params.M, params.N, params.fc/1e9, num_frames);
tag0 = sprintf('FAD2_N%d%s', params.N, ternary(quick, '_quick', ''));
diary([tag0 '_log.txt']); diary on; tic;
fprintf('%s\n', capt);

for s = 1:size(studies, 1)
    cfg = studies{s, 2}; nc = size(cfg, 1);
    bs = zeros(numel(snr_list), nc); bv = zeros(numel(vel_list), nc);
    for c = 1:nc
        pB = params; f = fieldnames(cfg{c,2});
        for q = 1:numel(f), pB.(f{q}) = cfg{c,2}.(f{q}); end
        pB.velocity_kmh = v_snr;
        t = sim_link(pB, snr_list, num_frames, 3000, false);  bs(:,c) = t(:,1);
        for i = 1:numel(vel_list)
            pB.velocity_kmh = vel_list(i);
            t = sim_link(pB, snr_for_vel, num_frames, 3000, false);  bv(i,c) = t(1);
        end
    end
    lbl = cfg(:,1); cols = lines(nc); st = arrayfun(@(k) {cols(k,:), '-o'}, 1:nc, 'UniformOutput', false);
    plain = regexprep(lbl, '[\\]', '');
    fprintf('\n[%s] BER vs SNR @ %d km/h\n SNR  %s\n', studies{s,1}, v_snr, sprintf('| %s ', plain{:}));
    for i = 1:numel(snr_list), fprintf('%4d  %s\n', snr_list(i), sprintf('%.2e  ', bs(i,:))); end
    fprintf('\n[%s] BER vs velocity @ %d dB\n km/h %s\n', studies{s,1}, snr_for_vel, sprintf('| %s ', plain{:}));
    for i = 1:numel(vel_list), fprintf('%4d  %s\n', vel_list(i), sprintf('%.2e  ', bv(i,:))); end

    fig = figure; h = plot_ber(snr_list(:), bs, nbits_total, 1e-5, st);
    xlabel('SNR (dB)'); ylabel('BER'); legend(h, lbl, 'Location', 'southwest', 'FontSize', 8);
    title({sprintf('Fading study %s: BER vs SNR (%d km/h)', studies{s,1}, v_snr), capt}, 'FontSize', 8);
    saveas(fig, sprintf('%s_%s_snr.png', tag0, studies{s,1}));
    fig = figure; h = plot_ber(vel_list(:), bv, nbits_total, 1e-5, st);
    xlabel('Velocity (km/h)'); ylabel('BER'); legend(h, lbl, 'Location', 'best', 'FontSize', 8);
    title({sprintf('Fading study %s: BER vs velocity (SNR %d dB)', studies{s,1}, snr_for_vel), capt}, 'FontSize', 8);
    saveas(fig, sprintf('%s_%s_vel.png', tag0, studies{s,1}));
    R.(studies{s,1}) = struct('labels', {plain}, 'snr_list', snr_list, 'bs', bs, 'vel_list', vel_list, 'bv', bv);
end
save([tag0 '_results.mat'], '-struct', 'R');
fprintf('\nTotal time: %.1f s\n', toc);
diary off;

function s = ternary(c, a, b)
    if c, s = a; else, s = b; end
end

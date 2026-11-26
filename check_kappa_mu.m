%% check_kappa_mu.m - (1) per-tap power/AF vs theory, (2) Monte-Carlo uncertainty of OTFS BER
% across kappa-mu(mu=1,2,4), Rician-all-taps and Rayleigh at 120 km/h, 15 dB.
clear; clc;
NFR = 600;                              % frames for part 2 (more = tighter error bars)
kappa = 10^(6/10);
cfg = {'rayleigh',        struct('fading_model','rayleigh');
       'rician-all-taps', struct('fading_model','rician-all-taps');
       'kappa-mu mu=1',   struct('fading_model','kappa-mu','kappa',kappa,'mu',1);
       'kappa-mu mu=2',   struct('fading_model','kappa-mu','kappa',kappa,'mu',2);
       'kappa-mu mu=4',   struct('fading_model','kappa-mu','kappa',kappa,'mu',4)};
mu_of = [NaN NaN 1 2 4];
base.M = 64; base.N = 8; base.delta_f = 15e3; base.fc = 4e9; base.velocity_kmh = 120; base.snr_db = 15;

%% 1. per-tap E|h|^2 and amount of fading (AF = var(|h|^2)/E^2)
fprintf('--- Part 1: per-tap statistics (20000 draws) ---\n');
for c = 1:size(cfg,1)
    p = base; f = fieldnames(cfg{c,2}); for q = 1:numel(f), p.(f{q}) = cfg{c,2}.(f{q}); end
    rng(3); G = zeros(20000, 9);
    for t = 1:20000, ch = channel_state(p); G(t,:) = abs(ch.h.').^2; end
    if c == 1, th = 1; elseif c == 2, th = (1+2*kappa)/(1+kappa)^2;
    else, th = (1+2*kappa)/(mu_of(c)*(1+kappa)^2); end
    fprintf('%-16s theory AF %.3f | measured AF min/max over taps %.3f / %.3f | E|h|^2 / p (tap1..9 ratio to Rayleigh-profile): %s\n', ...
        cfg{c,1}, th, min(var(G)./mean(G).^2), max(var(G)./mean(G).^2), mat2str(mean(G)./[0.2388 0.1717 0.1750 0.1055 0.2113 0.0295 0.0482 0.0151 0.0049], 3));
end

%% 2. BER with per-frame statistics (same seeds for every config)
fprintf('\n--- Part 2: OTFS BER, 120 km/h, 15 dB, %d frames ---\n', NFR);
nb = 2*base.M*base.N; fb = zeros(NFR, size(cfg,1));
for c = 1:size(cfg,1)
    p = base; f = fieldnames(cfg{c,2}); for q = 1:numel(f), p.(f{q}) = cfg{c,2}.(f{q}); end
    for fr = 1:NFR
        rng(3000 + fr); ch = channel_state(p); p.cp_len = ch.max_delay;
        bits = randi([0 1], nb, 1);
        tx = otfs_transmitter(p, bits);
        [y, H, nv] = otfs_channel(tx, p, ch);
        [~,~,~,e] = otfs_receiver(y, H, nv, bits, p);
        fb(fr, c) = e / nb;
    end
end
m = mean(fb); se = std(fb) / sqrt(NFR);
fprintf('%-16s BER        +/-1 SE     rel.SE   (95%% CI = +/-2 SE)\n', 'config');
for c = 1:size(cfg,1), fprintf('%-16s %.3e  %.2e  %5.1f%%\n', cfg{c,1}, m(c), se(c), 100*se(c)/m(c)); end
d = m(3) - m(2); sd = sqrt(se(3)^2 + se(2)^2);
fprintf('\nmu=1 minus Rician-all-taps (same distribution): diff %.2e, z = %.2f (|z|<2 => pure noise)\n', d, d/sd);
for c = 4:5
    d = m(c) - m(c-1); sd = sqrt(se(c)^2 + se(c-1)^2);
    fprintf('%s minus %s: diff %.2e, z = %.2f\n', cfg{c,1}, cfg{c-1,1}, d, d/sd);
end

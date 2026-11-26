%% crosscheck_toolbox.m - OPTIONAL: compare channel_state Rayleigh/Rician gain
% statistics with comm.RayleighChannel / comm.RicianChannel (Communications
% Toolbox). Statistical check only; the simulation keeps using channel_state.
% Compares, for tap 1 (largest-power EVA tap, normalised): mean power, amount
% of fading (var(|h|^2)/mean^2), and the Rician K-factor estimate.
clear; clc;
if isempty(ver('comm')), error('Communications Toolbox not found.'); end
rng(1);
T = 20000;                                   % samples / draws
params.M = 64; params.delta_f = 15e3; params.fc = 4e9;
params.velocity_kmh = 0;                     % Doppler irrelevant for gain stats
Kdb = 6; K = 10^(Kdb/10);
tau = [0 30 150 310 370 710 1090 1730 2510]; pdb = [0 -1.5 -1.4 -3.6 -0.6 -9.1 -7 -12 -16.9];
pp = 10.^(pdb/10); pp = pp/sum(pp);
% channel_state puts the whole LOS power K/(K+1) on tap 1, so tap-1 K is K/p1
Ktap1 = K / pp(1);
kest = @(g) kfac(abs(g).^2);

% --- our channel: tap 1, unit-power-normalised ---
mine = struct('rayleigh', zeros(T,1), 'rician', zeros(T,1));
for mdl = {'rayleigh','rician'}
    p = params; p.fading_model = mdl{1}; g = zeros(T,1);
    for t = 1:T, c = channel_state(p); g(t) = c.h(1); end
    mine.(mdl{1}) = g / sqrt(mean(abs(g).^2));
end

% --- toolbox: single-path, fixed (near-zero) Doppler, independent draws ---
% Large sample time gap decorrelates successive samples.
fs = 1e3; fd = 100;                          % 100 Hz Doppler, 1 kHz sampling
x = ones(T*50, 1);
ray = comm.RayleighChannel('SampleRate', fs, 'PathDelays', 0, ...
    'AveragePathGains', 0, 'MaximumDopplerShift', fd, 'PathGainsOutputPort', true);
[~, gR] = ray(x);  gR = gR(1:50:end);
ric = comm.RicianChannel('SampleRate', fs, 'PathDelays', 0, ...
    'AveragePathGains', 0, 'KFactor', Ktap1, 'MaximumDopplerShift', fd, ...
    'DirectPathDopplerShift', 0, 'PathGainsOutputPort', true);
[~, gC] = ric(x);  gC = gC(1:50:end);

fprintf('\nTap-1 gain statistics (Rayleigh AF=1; Rician K_tap1=%.2f -> AF=%.3f)\n', ...
    Ktap1, (1+2*Ktap1)/(1+Ktap1)^2);
fprintf('%-22s %-10s %-10s %-10s\n', 'source', 'E|h|^2', 'AF', 'K_est(dB)');
row = @(nm, g) fprintf('%-22s %-10.3f %-10.3f %-10.2f\n', nm, mean(abs(g).^2), ...
    var(abs(g).^2)/mean(abs(g).^2)^2, 10*log10(max(kest(g),eps)));
row('channel_state Rayleigh', mine.rayleigh);
row('comm.Rayleigh',          gR);
row('channel_state Rician',   mine.rician);
row('comm.Rician',            gC);
fprintf('Expect: E|h|^2 ~ 1, AF and K_est(dB) of each channel_state row within ~5%% of its comm.* row.\n');

function k = kfac(P)
% Moment-based K estimate from power samples: K = sqrt(1-AF')/(1-sqrt(1-AF')),
% AF' = var(P)/mean(P)^2 (valid for AF' < 1).
    a = var(P)/mean(P)^2;
    s = sqrt(max(1 - a, 0));
    k = s/(1 - s);
end

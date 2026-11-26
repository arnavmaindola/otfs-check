%% sanity_test.m - end-to-end checks of the OTFS/OFDM simulator (~1-2 min)
% Each check prints PASS/FAIL. Covers: transforms, power/SNR scaling,
% channel/matrix consistency, fading normalisation, MMSE correctness,
% data-leakage, and the expected BER trends.
clearvars -except N_GRID M_GRID; clc;
if ~exist('N_GRID','var'), N_GRID = 8; end    % set N_GRID=16 before running to test that grid
if ~exist('M_GRID','var'), M_GRID = 64; end   % set M_GRID=32 for the synopsis grid
params.M = M_GRID; params.N = N_GRID; params.delta_f = 15e3; params.fc = 4e9;
params.velocity_kmh = 120; params.snr_db = 10; params.fading_model = 'rayleigh';
M = params.M; N = params.N; MN = M*N;
nfail = 0;
pf = {'FAIL','PASS'};
chk = @(name, ok, info) fprintf('[%s] %-62s %s\n', pf{double(logical(ok))+1}, name, info);

%% 1. Transforms are unitary and FFT shortcut equals kron form
Fn = dftmtx(N)/sqrt(N); Fm = dftmtx(M)/sqrt(M);
ok = norm(Fn'*Fn - eye(N),'fro') < 1e-10 && norm(Fm'*Fm - eye(M),'fro') < 1e-10;
chk('DFT matrices unitary', ok, ''); nfail = nfail + ~ok;
rng(1); p = params; chan = channel_state(p);
[~, Hf] = otfs_channel(zeros(MN,1), p, chan);
U = kron(Fn, eye(M));
G = zeros(MN); n = (1:MN)';
for i = 1:numel(chan.h)
    idx = sub2ind([MN MN], n, mod(n-1-chan.lt(i), MN)+1);
    G(idx) = G(idx) + chan.h(i)*exp(1i*2*pi*chan.nu(i)*(n-1)*chan.Ts);
end
e = norm(Hf - U*G*U','fro')/norm(Hf,'fro');
chk('H_eff (FFT form) == U*G*U'' (kron form)', e < 1e-10, sprintf('err=%.1e', e)); nfail = nfail + ~(e<1e-10);

%% 2. Tx power, loop-back with ideal channel (no leakage path needed)
bits = randi([0 1], 2*MN, 1); p.cp_len = 3;
[tx, X] = otfs_transmitter(p, bits); [txo, Xo] = ofdm_transmitter(p, bits);
pw1 = mean(abs(tx).^2); pw2 = mean(abs(txo).^2);
ok = abs(pw1-1) < 0.1 && abs(pw2-1) < 0.1;
chk('Unit average Tx power (OTFS, OFDM)', ok, sprintf('%.3f, %.3f', pw1, pw2)); nfail = nfail + ~ok;
ok = isequal(X, Xo);
chk('OTFS and OFDM carry identical symbol grids', ok, ''); nfail = nfail + ~ok;
Xr = fft(reshape(tx, M, N), [], 2)/sqrt(N);
ok = norm(Xr - X,'fro') < 1e-10;
chk('OTFS Tx -> SFFT recovers DD grid exactly', ok, ''); nfail = nfail + ~ok;
[rb, ~, ber] = otfs_receiver(X, eye(MN), 1e-6, bits, p);
ok = ber == 0 && isequal(rb, bits);
chk('OTFS loop-back (H=I, no noise) BER = 0', ok, sprintf('BER=%g', ber)); nfail = nfail + ~ok;
[~, ~, ber] = ofdm_receiver(Xo, eye(MN), 1e-6, bits, p, 'mmse');
chk('OFDM loop-back (H=I, no noise) BER = 0', ber==0, sprintf('BER=%g', ber)); nfail = nfail + ~(ber==0);

%% 3. Noise variance / SNR definition
p2 = params; p2.snr_db = 10; rng(2);
[y, ~, nv] = otfs_channel(zeros(MN,1), p2, chan);
ok = abs(mean(abs(y(:)).^2)/nv - 1) < 0.2;
chk('OTFS noise power equals 10^(-SNR/10)', ok, sprintf('ratio=%.2f', mean(abs(y(:)).^2)/nv)); nfail = nfail + ~ok;
p2.cp_len = chan.max_delay;
[y, ~, nv] = ofdm_channel(zeros((M+p2.cp_len)*N,1), p2, chan);
ok = abs(mean(abs(y(:)).^2)/nv - 1) < 0.2;
chk('OFDM noise power equals 10^(-SNR/10)', ok, sprintf('ratio=%.2f', mean(abs(y(:)).^2)/nv)); nfail = nfail + ~ok;

%% 4. Noise-free channel output equals H * x (both systems)
p3 = params; p3.snr_db = 300; p3.cp_len = chan.max_delay;
[txa, Xa] = otfs_transmitter(p3, bits);
[y, H] = otfs_channel(txa, p3, chan);
e = norm(y(:) - H*Xa(:))/norm(y(:));
chk('OTFS: rx = H_eff * x (noise-free)', e < 1e-8, sprintf('err=%.1e', e)); nfail = nfail + ~(e<1e-8);
[txb, Xb] = ofdm_transmitter(p3, bits);
[y, H] = ofdm_channel(txb, p3, chan);
e = norm(y(:) - H*Xb(:))/norm(y(:));
chk('OFDM: rx = H_tf * x (noise-free)', e < 1e-8, sprintf('err=%.1e', e)); nfail = nfail + ~(e<1e-8);

%% 5. OFDM ICI: ~0 with CP and no Doppler, grows with velocity
ici = zeros(1,2);
for k = 1:2
    vv = [0 500]; pv = params; pv.velocity_kmh = vv(k); rng(5);
    c = channel_state(pv); pv.cp_len = c.max_delay;
    pv.snr_db = 300; [~, Ht] = ofdm_channel(zeros((M+c.max_delay)*N,1), pv, c);
    ici(k) = 1 - sum(abs(diag(Ht)).^2)/sum(abs(Ht(:)).^2);
end
ok = ici(1) < 1e-10 && ici(2) > 1e-3;
chk('OFDM ICI fraction: 0 km/h ~ 0, 500 km/h > 0', ok, sprintf('%.1e, %.1e', ici)); nfail = nfail + ~ok;

%% 6. Fading models: E[sum|h|^2] = 1 and amount-of-fading ordering
models = {'rayleigh','rician','nakagami-m','alpha-mu','kappa-mu'};
pwr = zeros(1,5); AF = zeros(1,5); rng(7);
for m = 1:5
    pm = params; pm.fading_model = models{m}; g = zeros(3000,1);
    for t = 1:3000, c = channel_state(pm); g(t) = abs(c.h(1))^2; pwr(m) = pwr(m) + sum(abs(c.h).^2)/3000; end
    AF(m) = var(g)/mean(g)^2;   % amount of fading of the strongest tap
end
ok = all(abs(pwr-1) < 0.08);
chk('All 5 models power-normalised (E sum|h|^2 = 1)', ok, sprintf('%.2f ', pwr)); nfail = nfail + ~ok;
ok = AF(4) > AF(1) && AF(1) > AF(3) && AF(1) > AF(2) && AF(1) > AF(5);
chk('Fading severity: a-mu > Rayleigh > Rician,Nakagami,k-mu', ok, sprintf('AF=%.2f ', AF)); nfail = nfail + ~ok;

%% 7. MMSE correctness and no data leakage
rng(11); p = params; p.snr_db = 5; c = channel_state(p); p.cp_len = c.max_delay;
bits = randi([0 1], 2*MN, 1);
[tx, Xg] = otfs_transmitter(p, bits); [y, H, nv] = otfs_channel(tx, p, c);
[rb1, xs] = otfs_receiver(y, H, nv, bits, p);
[rb2] = otfs_receiver(y, H, nv, 1-bits, p);          % wrong "ground truth"
ok = isequal(rb1, rb2);
chk('Receiver output independent of tx_bits (no leakage)', ok, ''); nfail = nfail + ~ok;
xm = H' * ((H*H' + nv*eye(MN)) \ y(:));              % alternative MMSE form
e = norm(xm - xs(:))/norm(xm);
chk('MMSE: (H''H+s2 I)^-1 H'' y == H''(HH''+s2 I)^-1 y', e < 1e-6, sprintf('err=%.1e', e)); nfail = nfail + ~(e<1e-6);
% MMSE <= ZF in MSE at low SNR
mse_m = 0; mse_z = 0;
for t = 1:15
    rng(100+t); c = channel_state(p); [tx, Xg] = otfs_transmitter(p); [y, H, nv] = otfs_channel(tx, p, c);
    [~, xs] = otfs_receiver(y, H, nv, zeros(2*MN,1), p);
    xz = H \ y(:); mse_m = mse_m + norm(xs(:)-Xg(:))^2; mse_z = mse_z + norm(xz-Xg(:))^2;
end
ok = mse_m < mse_z;
chk('MMSE MSE < ZF MSE at 5 dB', ok, sprintf('%.1f < %.1f', mse_m, mse_z)); nfail = nfail + ~ok;
% pure-noise run must give BER ~ 0.5 (a leak would push it below)
pn = p; pn.snr_db = -40; b = sim_link(pn, -40, 5, 50, true);
ok = all(abs(b - 0.5) < 0.05);
chk('BER ~ 0.5 at -40 dB (no information leaking to Rx)', ok, sprintf('%.3f ', b)); nfail = nfail + ~ok;

%% 8. Trends (fast, 40 frames, common random numbers)
nf = 40; pt = params;
pt.velocity_kmh = 500;
b = sim_link(pt, [0 10 20 30], nf, 700, true);
ok = all(diff(b(:,1)) <= 0) && b(1,1) > 100*b(4,1) && all(diff(b(:,2)) < 0);
chk('BER falls with SNR (OTFS and OFDM, 500 km/h)', ok, ''); nfail = nfail + ~ok;
ok = b(3,1) < b(3,2) && b(4,1) < b(4,2);
chk('OTFS beats OFDM one-tap at 20/30 dB, 500 km/h', ok, sprintf('%.1e vs %.1e @30dB', b(4,1), b(4,2))); nfail = nfail + ~ok;
ok = b(4,2) > 3*b(4,3);
chk('OFDM one-tap shows ICI error floor vs full MMSE @30 dB', ok, sprintf('%.1e vs %.1e', b(4,2), b(4,3))); nfail = nfail + ~ok;
bo = zeros(2,3);
vv2 = [0 700];
for k = 1:2, pt.velocity_kmh = vv2(k); bo(k,:) = sim_link(pt, 25, nf, 800, true); end
ok = bo(2,2) > 2*bo(1,2);
chk('OFDM BER grows with velocity (0 -> 700 km/h, 25 dB)', ok, sprintf('%.1e -> %.1e', bo(1,2), bo(2,2))); nfail = nfail + ~ok;
ok = bo(2,1) < 2*bo(1,1) + 1e-3;
chk('OTFS BER does not degrade with velocity (25 dB)', ok, sprintf('%.1e -> %.1e', bo(1,1), bo(2,1))); nfail = nfail + ~ok;

fprintf('\n%d check(s) failed.\n', nfail);


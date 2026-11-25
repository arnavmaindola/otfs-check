%% main2.m - OTFS vs OFDM, and OTFS under 5 fading models
% Part A : OTFS vs OFDM   (BER vs SNR, BER vs velocity)
% Part B : OTFS in Rayleigh / Rician / Nakagami-m / alpha-mu / kappa-mu
% Run time ~ 5-8 min at num_frames = 100 (set quick = true for ~1 min).
clear; clc; close all;

quick = false;
num_frames = 100; if quick, num_frames = 20; end

params.M = 64;  params.N = 8;          % MN = 512 (delay x Doppler grid)
params.delta_f = 15e3;  params.fc = 4e9;
params.velocity_kmh = 120; params.snr_db = 10; params.fading_model = 'rayleigh';

snr_list = 0:5:30;
vel_list = 0:100:700;
snr_for_vel = 20;                      % Part A  BER-vs-velocity SNR
snr_for_vel_B = 15;                    % Part B  BER-vs-velocity SNR
v_B = 120;                             % Part B  BER-vs-SNR velocity
models = {'rayleigh','rician','nakagami-m','alpha-mu','kappa-mu'};
names  = {'Rayleigh','Rician (K=6dB)','Nakagami-m (m=2)','\alpha-\mu (1.2,0.8)','\kappa-\mu (3,2)'};
mk = {'o-','s-','^-','d-','v-'};
floor_ber = 0.5 / (2*params.M*params.N*num_frames);   % plot floor for 0 errors
pl = @(b) max(b, floor_ber);

%% ---------------- Part A : OTFS vs OFDM ----------------
tic;
figure('Position',[50 50 1100 450]);
for k = 1:2
    v = [120 500]; pA = params; pA.velocity_kmh = v(k);
    b = sim_link(pA, snr_list, num_frames, 1000, true);
    fprintf('\nA1: BER vs SNR @ %d km/h\n SNR   OTFS     OFDM-1tap OFDM-MMSE\n', v(k));
    fprintf('%4d  %.2e  %.2e  %.2e\n', [snr_list(:) b]');
    subplot(1,2,k);
    semilogy(snr_list, pl(b(:,1)), 'b-o', snr_list, pl(b(:,2)), 'r-s', ...
             snr_list, pl(b(:,3)), 'm--^', 'LineWidth', 1.5); grid on;
    xlabel('SNR (dB)'); ylabel('BER'); ylim([floor_ber*0.8 1]);
    title(sprintf('OTFS vs OFDM, %d km/h', v(k)));
    legend('OTFS (MMSE)','OFDM (one-tap)','OFDM (full MMSE)','Location','southwest');
end

bv = zeros(numel(vel_list), 3);
for i = 1:numel(vel_list)
    pA = params; pA.velocity_kmh = vel_list(i);
    bv(i,:) = sim_link(pA, snr_for_vel, num_frames, 1000, true);
end
fprintf('\nA2: BER vs velocity @ %d dB\n km/h   OTFS     OFDM-1tap OFDM-MMSE\n', snr_for_vel);
fprintf('%4d  %.2e  %.2e  %.2e\n', [vel_list(:) bv]');
figure;
semilogy(vel_list, pl(bv(:,1)), 'b-o', vel_list, pl(bv(:,2)), 'r-s', ...
         vel_list, pl(bv(:,3)), 'm--^', 'LineWidth', 1.5); grid on;
xlabel('Velocity (km/h)'); ylabel('BER'); ylim([floor_ber*0.8 1]);
title(sprintf('OTFS vs OFDM: BER vs velocity (SNR = %d dB)', snr_for_vel));
legend('OTFS (MMSE)','OFDM (one-tap)','OFDM (full MMSE)','Location','northwest');

%% ---------------- Part B : fading models (OTFS only) ----------------
bs = zeros(numel(snr_list), 5);  bvB = zeros(numel(vel_list), 5);
for m = 1:5
    pB = params; pB.fading_model = models{m};
    pB.velocity_kmh = v_B;
    bs(:,m) = sim_link(pB, snr_list, num_frames, 3000, false);
    for i = 1:numel(vel_list)
        pB.velocity_kmh = vel_list(i);
        bvB(i,m) = sim_link(pB, snr_for_vel_B, num_frames, 3000, false);
    end
end
fprintf('\nB1: OTFS BER vs SNR @ %d km/h\n SNR  %s\n', v_B, sprintf('%-11s', models{:}));
for i = 1:numel(snr_list), fprintf('%4d  %s\n', snr_list(i), sprintf('%.2e   ', bs(i,:))); end
fprintf('\nB2: OTFS BER vs velocity @ %d dB\n km/h %s\n', snr_for_vel_B, sprintf('%-11s', models{:}));
for i = 1:numel(vel_list), fprintf('%4d  %s\n', vel_list(i), sprintf('%.2e   ', bvB(i,:))); end

figure; hold on;
for m = 1:5, semilogy(snr_list, pl(bs(:,m)), mk{m}, 'LineWidth', 1.5); end
set(gca,'YScale','log'); grid on; ylim([floor_ber*0.8 1]);
xlabel('SNR (dB)'); ylabel('BER'); legend(names, 'Location','southwest');
title(sprintf('OTFS BER vs SNR, fading models (%d km/h)', v_B));

figure; hold on;
for m = 1:5, semilogy(vel_list, pl(bvB(:,m)), mk{m}, 'LineWidth', 1.5); end
set(gca,'YScale','log'); grid on; ylim([floor_ber*0.8 1]);
xlabel('Velocity (km/h)'); ylabel('BER'); legend(names, 'Location','best');
title(sprintf('OTFS BER vs velocity, fading models (SNR = %d dB)', snr_for_vel_B));
fprintf('\nTotal time: %.1f s\n', toc);

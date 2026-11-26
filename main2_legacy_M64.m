%% main2.m - OTFS vs OFDM, and OTFS under 5 fading models
% Part A : OTFS vs OFDM   (BER vs SNR, BER vs velocity)
% Part B : OTFS in Rayleigh / Rician / Nakagami-m / alpha-mu / kappa-mu
% Run time ~ 5-8 min at num_frames = 100 (set quick = true for ~1 min).
clear; clc; close all;

quick = false;
NUM_FRAMES = 100;                      % Monte-Carlo frames per point; use 300 for the final run (~4.5 min)
num_frames = NUM_FRAMES; if quick, num_frames = 20; end

N_GRID = 8;                            % Doppler bins N (grid is M x N; 16 -> see main_n16)
params.M = 64;  params.N = N_GRID;     % MN = 512 for N = 8
params.delta_f = 15e3;  params.fc = 4e9;
params.velocity_kmh = 120; params.snr_db = 10; params.fading_model = 'rayleigh';

snr_list = 0:5:30;
vel_list = 0:50:700;                   % Part A velocity grid (Part B uses vel_list_B)
vel_list_B = 0:100:700;
v_list_A = [120 300 500];              % Part A BER-vs-SNR velocities
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
run_partA(params, num_frames, snr_list, v_list_A, vel_list, snr_for_vel, sprintf('A_N%d', N_GRID));

%% ---------------- Part B : fading models (OTFS only) ----------------
bs = zeros(numel(snr_list), 5);  bvB = zeros(numel(vel_list_B), 5);
for m = 1:5
    pB = params; pB.fading_model = models{m};
    pB.velocity_kmh = v_B;
    t = sim_link(pB, snr_list, num_frames, 3000, false);  bs(:,m) = t(:,1);
    for i = 1:numel(vel_list_B)
        pB.velocity_kmh = vel_list_B(i);
        t = sim_link(pB, snr_for_vel_B, num_frames, 3000, false);  bvB(i,m) = t(1);
    end
end
fprintf('\nB1: OTFS BER vs SNR @ %d km/h\n SNR  %s\n', v_B, sprintf('%-11s', models{:}));
for i = 1:numel(snr_list), fprintf('%4d  %s\n', snr_list(i), sprintf('%.2e   ', bs(i,:))); end
fprintf('\nB2: OTFS BER vs velocity @ %d dB\n km/h %s\n', snr_for_vel_B, sprintf('%-11s', models{:}));
for i = 1:numel(vel_list_B), fprintf('%4d  %s\n', vel_list_B(i), sprintf('%.2e   ', bvB(i,:))); end

figure; hold on;
for m = 1:5, semilogy(snr_list, pl(bs(:,m)), mk{m}, 'LineWidth', 1.5); end
set(gca,'YScale','log'); grid on; ylim([floor_ber*0.8 1]);
xlabel('SNR (dB)'); ylabel('BER'); legend(names, 'Location','southwest');
title(sprintf('OTFS BER vs SNR, fading models (%d km/h)', v_B));

figure; hold on;
for m = 1:5, semilogy(vel_list_B, pl(bvB(:,m)), mk{m}, 'LineWidth', 1.5); end
set(gca,'YScale','log'); grid on; ylim([floor_ber*0.8 1]);
xlabel('Velocity (km/h)'); ylabel('BER'); legend(names, 'Location','best');
title(sprintf('OTFS BER vs velocity, fading models (SNR = %d dB)', snr_for_vel_B));
fprintf('\nTotal time: %.1f s\n', toc);

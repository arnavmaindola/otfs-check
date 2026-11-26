%% main2.m - FINAL Part A: OTFS vs OFDM over EVA (M=32, N=16 per synopsis)
% BER vs SNR at 120/300/500 km/h and BER vs velocity (0:50:500 km/h, 20 dB).
% Part B now lives in main_partB_final.m. Outputs are tagged A_M32_N16*.
% (Previous M=64 runner kept as main2_legacy_M64.m.)
% Runtime estimate at 1000 frames: ~40-60 min; run quick = true first (~2 min).
clear; clc; close all;

quick = true;                          % true: 20 frames, coarse grids (smoke test)
NUM_FRAMES = 1000;                     % Monte-Carlo frames for final results
N_GRID = 16;  M_GRID = 32;             % OTFS/OFDM grid (synopsis: M=32, N=16)
num_frames = NUM_FRAMES; if quick, num_frames = 20; end

params.M = M_GRID;  params.N = N_GRID;
params.delta_f = 15e3;  params.fc = 4e9;
params.velocity_kmh = 120; params.snr_db = 10; params.fading_model = 'rayleigh';

snr_list = 0:5:30;
vel_list = 0:50:500;
v_list_A = [120 300 500];              % BER-vs-SNR velocities
snr_for_vel = 20;                      % BER-vs-velocity SNR
if quick, snr_list = 0:10:30; vel_list = [0 250 500]; end

tic;
run_partA(params, num_frames, snr_list, v_list_A, vel_list, snr_for_vel, ...
          sprintf('A_M%d_N%d%s', M_GRID, N_GRID, ternary(quick, '_quick', '')));
fprintf('\nTotal time: %.1f s\n', toc);

function s = ternary(c, a, b)
    if c, s = a; else, s = b; end
end

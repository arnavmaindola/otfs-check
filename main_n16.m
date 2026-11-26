%% main_n16.m - OTFS vs OFDM on a 64 x 16 grid (Part A only, reduced cost)
% N = 16 -> MN = 1024: every MMSE solve is ~8x costlier than at N = 8.
% Reduced sweep: 100 frames, SNR {0,10,20,30} dB, velocities {120,300,500},
% BER-vs-velocity at 20 dB over 0:100:700. Expect roughly 5-15 min on a
% desktop and longer on MATLAB Online; set quick = true (20 frames) first.
% Outputs are prefixed A_N16_ so the N = 8 run (A_N8_*) is not overwritten.
clear; clc; close all;

quick = false;
NUM_FRAMES = 100;                      % Monte-Carlo frames per point
num_frames = NUM_FRAMES; if quick, num_frames = 20; end
N_GRID = 16;

params.M = 64;  params.N = N_GRID;
params.delta_f = 15e3;  params.fc = 4e9;
params.velocity_kmh = 120; params.snr_db = 10; params.fading_model = 'rayleigh';

snr_list = 0:10:30;
v_list_A = [120 300 500];
vel_list = 0:100:700;
snr_for_vel = 20;

tic;
run_partA(params, num_frames, snr_list, v_list_A, vel_list, snr_for_vel, sprintf('A_N%d', N_GRID));
fprintf('\nTotal time: %.1f s\n', toc);

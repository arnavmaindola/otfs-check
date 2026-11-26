%% plot_dd_channel.m - delay-Doppler response of the OTFS effective channel
% Shows |H_eff| response to a unit impulse at DD (l=0,k=0) for ONE channel
% draw (same seed => same fades and path angles) at several speeds. Doppler
% spread grows with speed. Dots mark the true (delay, Doppler) of each path.
clear; clc; close all;

N_GRID = 16;  M_GRID = 32;             % synopsis grid
vel_list = [0 300 500];
seed = 1000 + 1;
l_show = 0:3;                          % delay taps shown (EVA max delay is 1 tap at M=32)

params.M = M_GRID;  params.N = N_GRID;  params.delta_f = 15e3;  params.fc = 4e9;
params.snr_db = 30; params.fading_model = 'rayleigh';
M = params.M; N = params.N; Ts = 1/(M*params.delta_f);
dop_bins = -N/2 : N/2-1;                          % centred Doppler bins
dop_hz = dop_bins * params.delta_f / N;           % bin spacing delta_f/N

fig = figure('Position', [30 60 430*numel(vel_list) 420]);
for iv = 1:numel(vel_list)
    p = params; p.velocity_kmh = vel_list(iv);
    rng(seed); chan = channel_state(p);
    [~, H] = otfs_channel(zeros(M*N,1), p, chan);
    R = abs(reshape(H(:,1), M, N));               % response at (delay l, Doppler k)
    R = fftshift(R, 2);                           % centre Doppler
    subplot(1, numel(vel_list), iv);
    imagesc(dop_hz/1e3, l_show*Ts*1e6, R(l_show+1, :)); axis xy; colorbar; hold on;
    plot(chan.nu/1e3, chan.lt*Ts*1e6, 'wo', 'MarkerFaceColor', 'r');
    xlim([dop_hz(1) - p.delta_f/N/2, dop_hz(end) + p.delta_f/N/2]/1e3);
    xlabel('Doppler (kHz); bin width \Deltaf/N'); ylabel('Delay (\mus)');
    title({sprintf('%d km/h (max Doppler %.0f Hz)', vel_list(iv), ...
           (vel_list(iv)/3.6)*params.fc/3e8), ...
           sprintf('M=%d, N=%d, f_c=%g GHz', M, N, params.fc/1e9)}, 'FontSize', 9);
end
sgtitle('OTFS effective DD channel |H_{eff}| (impulse response, same draw across speeds); dots = path delay/Doppler');
saveas(fig, sprintf('DD_channel_M%d_N%d.png', M, N));

function res = run_partA(params, num_frames, snr_list, v_snr_list, vel_list, snr_for_vel, tag)
% RUN_PARTA  OTFS vs OFDM: BER vs SNR (one panel per velocity in v_snr_list)
% and BER vs velocity. Prints tables, draws figures, and saves everything
% under the prefix 'tag' (<tag>_log.txt, <tag>_results.mat, <tag>_*.png).
%   Dashed black = analytic OFDM one-tap approximation (invalid at <=5 dB: ignores noise-limited regime)
%   BER ~ 0.5*(sigma_ici^2 + 1/SNR), sigma_ici^2 = (pi*fd/delta_f)^2/6.
%   Hollow markers = zero errors observed (BER < 1/bits, drawn at the floor).

    M = params.M; N = params.N;
    ylo = 5e-7;   % axis bottom; real floor is 1/bits (~9.8e-7 at 1000 frames, M=32, N=16)
    nbits_total = 2 * M * N * num_frames;
    capt = sprintf('M=%d N=%d fc=%gGHz %d frames', M, N, params.fc/1e9, num_frames);
    st = {{'b','-o'}, {'r','-s'}, {'m','--^'}};
    lg = {'OTFS (MMSE)','OFDM (one-tap)','OFDM (full MMSE)'};

    diary([tag '_log.txt']); diary on;
    fprintf('%s  (hollow = 0 errors, floor 1/%d = %.1e)\n', capt, nbits_total, 1/nbits_total);

    nv = numel(v_snr_list);
    bsnr = zeros(numel(snr_list), 3, nv);
    fig1 = figure('Position', [30 50 480*nv 480]);
    for k = 1:nv
        pA = params; pA.velocity_kmh = v_snr_list(k);
        b = sim_link(pA, snr_list, num_frames, 1000, true);
        bsnr(:,:,k) = b;
        fprintf('\nA1: BER vs SNR @ %d km/h\n SNR   OTFS     OFDM-1tap OFDM-MMSE  OFDM-1tap(analytic)\n', v_snr_list(k));
        fd = v_snr_list(k)/3.6 * params.fc / 3e8;
        s_ici = (pi*fd/params.delta_f)^2 / 6;
        ban = min(0.5*(s_ici + 10.^(-snr_list(:)/10)), 0.5);
        fprintf('%4d  %.2e  %.2e  %.2e  %.2e\n', [snr_list(:) b ban]');
        subplot(1, nv, k);
        h = plot_ber(snr_list(:), b, nbits_total, ylo, st);
        ha = semilogy(snr_list, ban, 'k--', 'LineWidth', 1.2);
        xlabel('SNR (dB)'); ylabel('BER');
        title({sprintf('OTFS vs OFDM, %d km/h', v_snr_list(k)), capt}, 'FontSize', 8);
        legend([h; ha], [lg, {'OFDM 1-tap analytic (valid >5 dB)'}], 'Location', 'southwest', 'FontSize', 8);
    end
    saveas(fig1, [tag '_A1_ber_vs_snr.png']);

    bv = zeros(numel(vel_list), 3);
    for i = 1:numel(vel_list)
        pA = params; pA.velocity_kmh = vel_list(i);
        bv(i,:) = sim_link(pA, snr_for_vel, num_frames, 1000, true);
    end
    fprintf('\nA2: BER vs velocity @ %d dB\n km/h   OTFS     OFDM-1tap OFDM-MMSE\n', snr_for_vel);
    fprintf('%4d  %.2e  %.2e  %.2e\n', [vel_list(:) bv]');
    fig2 = figure;
    h = plot_ber(vel_list(:), bv, nbits_total, ylo, st);
    xlabel('Velocity (km/h)'); ylabel('BER');
    title({sprintf('OTFS vs OFDM: BER vs velocity (SNR = %d dB)', snr_for_vel), capt}, 'FontSize', 8);
    legend(h, lg, 'Location', 'northwest');
    saveas(fig2, [tag '_A2_ber_vs_velocity.png']);

    res = struct('snr_list', snr_list, 'v_snr_list', v_snr_list, 'bsnr', bsnr, ...
                 'vel_list', vel_list, 'snr_for_vel', snr_for_vel, 'bv', bv, ...
                 'num_frames', num_frames, 'M', M, 'N', N);
    save([tag '_results.mat'], '-struct', 'res');
    diary off;
end



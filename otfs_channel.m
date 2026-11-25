function [rx_grid_dd, H_eff, noise_var] = otfs_channel(tx_sig, params, chan)
% OTFS channel: EVA multipath + Doppler + fading + AWGN.
%
% Inputs:
%   tx_sig - Transmitted time-domain signal (MN x 1)
%   params - M, N, delta_f, fc, velocity_kmh, snr_db, fading_model
%   chan   - (optional) channel realisation from channel_state(); drawn
%            internally if omitted. Pass the same 'chan' to ofdm_channel
%            for a fair OTFS/OFDM comparison.
%
% Outputs:
%   rx_grid_dd - Received delay-Doppler grid (M x N)
%   H_eff      - Effective DD-domain channel (MN x MN) = U*G*U'
%   noise_var  - Noise variance (SNR = 1/noise_var, unit-power signal)
%
% One frame-level cyclic prefix is assumed, so the time-domain channel G
% is circulant over the MN samples of the frame (the usual OTFS model).

    M = params.M; N = params.N; MN = M*N;
    if nargin < 3 || isempty(chan), chan = channel_state(params); end

    G = zeros(MN, MN);
    n = (1:MN).';
    for i = 1:numel(chan.h)
        cols = mod(n - 1 - chan.lt(i), MN) + 1;
        idx  = sub2ind([MN MN], n, cols);
        G(idx) = G(idx) + chan.h(i) * exp(1i*2*pi*chan.nu(i)*(n-1)*chan.Ts);
    end

    % H_eff = U*G*U', U = kron(F_N, I_M) (DFT along Doppler), via FFTs
    G4 = reshape(G, M, N, M, N);
    H_eff = reshape(ifft(fft(G4, [], 2), [], 4), MN, MN);

    noise_var = 10^(-params.snr_db / 10);

    y_t = reshape(G * tx_sig(:), M, N);
    y_dd = fft(y_t, [], 2) / sqrt(N);
    w = sqrt(noise_var/2) * (randn(M, N) + 1i*randn(M, N));
    rx_grid_dd = y_dd + w;
end

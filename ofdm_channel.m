function [rx_grid_tf, H_tf, noise_var] = ofdm_channel(tx_sig, params, chan)
% OFDM channel: same EVA/Doppler/fading realisation as otfs_channel.
% The CP'd signal passes through the (linear) time-varying channel, the CP
% is discarded, AWGN is added and an FFT gives the time-frequency grid.
%
% Outputs:
%   rx_grid_tf - Received TF grid (M x N)
%   H_tf       - Full TF-domain channel (MN x MN). Its diagonal is the
%                per-subcarrier gain; off-diagonals are the Doppler ICI.
%   noise_var  - Noise variance per resource element

    M = params.M; N = params.N; Lcp = params.cp_len;
    L = M + Lcp; Lt = N * L; MN = M * N;
    if nargin < 3 || isempty(chan), chan = channel_state(params); end
    if chan.max_delay > Lcp
        error('CP (%d) shorter than channel delay spread (%d)', Lcp, chan.max_delay);
    end

    Gl = zeros(Lt, Lt);                         % linear time-varying channel
    for i = 1:numel(chan.h)
        r = (chan.lt(i)+1 : Lt).';
        idx = sub2ind([Lt Lt], r, r - chan.lt(i));
        Gl(idx) = Gl(idx) + chan.h(i) * exp(1i*2*pi*chan.nu(i)*(r-1)*chan.Ts);
    end

    % Body-sample rows, and CP insertion matrix T (MN -> Lt)
    keep = reshape((1:Lt).', L, N);  keep = keep(Lcp+1:end, :);  keep = keep(:);
    T = sparse(Lt, MN);
    for k = 1:N
        body = (k-1)*L + Lcp + (1:M);
        T(body, (k-1)*M + (1:M)) = speye(M);
        if Lcp > 0
            T((k-1)*L + (1:Lcp), (k-1)*M + (M-Lcp+1:M)) = speye(Lcp);
        end
    end
    A = reshape(Gl(keep, :) * T, M, N, M, N);
    H_tf = reshape(ifft(fft(A, [], 1), [], 3), MN, MN);   % (I_N kron F_M) A (I_N kron F_M')

    noise_var = 10^(-params.snr_db / 10);

    y_t = reshape(Gl * tx_sig(:), L, N);
    y_t = y_t(Lcp+1:end, :);
    w = sqrt(noise_var/2) * (randn(M, N) + 1i*randn(M, N));
    rx_grid_tf = fft(y_t, [], 1) / sqrt(M) + w;   % unitary FFT keeps noise white
end

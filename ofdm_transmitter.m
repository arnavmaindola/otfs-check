function [tx_sig, tx_grid_tf, tx_bits] = ofdm_transmitter(params, tx_bits)
% OFDM transmitter: M subcarriers x N symbols, per-symbol cyclic prefix.
% Same QPSK mapping/grid as otfs_transmitter, so identical bits can be
% sent through both systems.
%
% params.cp_len : CP length in samples (>= channel max delay)
% tx_sig        : ((M+cp_len)*N x 1), unit average power per sample

    M = params.M; N = params.N; Lcp = params.cp_len;

    if nargin < 2 || isempty(tx_bits)
        tx_bits = randi([0, 1], M * N * 2, 1);
    end

    b0 = tx_bits(1:2:end);
    b1 = tx_bits(2:2:end);
    tx_syms = ((1 - 2*b0) + 1i*(1 - 2*b1)) / sqrt(2);
    tx_grid_tf = reshape(tx_syms, M, N);

    s = ifft(tx_grid_tf, [], 1) * sqrt(M);        % unitary IFFT per symbol
    s = [s(M-Lcp+1:M, :); s];                     % add CP
    tx_sig = s(:);
end

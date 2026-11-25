function [tx_sig, tx_grid_dd, tx_bits] = otfs_transmitter(params, tx_bits)
% OTFS Transmitter
% Generates Delay-Doppler QPSK symbols and modulates to time domain.
%
% Inputs:
%   params  - Simulation parameters struct (M, N)
%   tx_bits - Binary data bits (optional, generated randomly if empty)
%
% Outputs:
%   tx_sig     - Transmitted time-domain signal (MN x 1)
%   tx_grid_dd - Delay-Doppler symbol grid (M x N)
%   tx_bits    - Transmitted binary bits

    M = params.M;
    N = params.N;

    % Generate random bits if not provided (2 bits per QPSK symbol)
    if nargin < 2 || isempty(tx_bits)
        tx_bits = randi([0, 1], M * N * 2, 1);
    end

    % QPSK modulation with Gray mapping and unit average energy
    b0 = tx_bits(1:2:end);
    b1 = tx_bits(2:2:end);
    tx_syms = ((1 - 2*b0) + 1i*(1 - 2*b1)) / sqrt(2);

    % Map symbols into M x N Delay-Doppler grid
    tx_grid_dd = reshape(tx_syms, M, N);

    % Modulate to time domain: ISFFT + Heisenberg transform
    FN = dftmtx(N) / sqrt(N);
    S_mat = tx_grid_dd * FN';
    tx_sig = S_mat(:);
end

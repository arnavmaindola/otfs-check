function [rx_bits, rx_syms_dd, ber, num_errors] = otfs_receiver(rx_grid_dd, H_eff, noise_var, tx_bits, params)
% OTFS Receiver
% Performs MMSE equalization in Delay-Doppler domain and QPSK demodulation.
%
% Inputs:
%   rx_grid_dd - Received Delay-Doppler grid (M x N)
%   H_eff      - Effective channel matrix (MN x MN)
%   noise_var  - Noise variance sigma^2
%   tx_bits    - Transmitted bits (for BER computation)
%   params     - Simulation parameters struct (M, N)
%
% Outputs:
%   rx_bits    - Detected binary bits
%   rx_syms_dd - Equalized Delay-Doppler symbols (M x N)
%   ber        - Bit error rate
%   num_errors - Total bit errors

    M = params.M;
    N = params.N;
    MN = M * N;

    y_dd = rx_grid_dd(:);

    % MMSE equalization: x_hat = (H'*H + sigma^2*I) \ (H'*y)
    x_hat = (H_eff' * H_eff + noise_var * eye(MN)) \ (H_eff' * y_dd);
    rx_syms_dd = reshape(x_hat, M, N);

    % QPSK Gray demapping
    b0_hat = double(real(x_hat) < 0);
    b1_hat = double(imag(x_hat) < 0);
    rx_bits = zeros(2 * MN, 1);
    rx_bits(1:2:end) = b0_hat;
    rx_bits(2:2:end) = b1_hat;

    % BER calculation
    num_errors = sum(tx_bits(:) ~= rx_bits(:));
    ber = num_errors / length(tx_bits);
end

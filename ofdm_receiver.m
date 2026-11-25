function [rx_bits, rx_syms, ber, num_errors] = ofdm_receiver(rx_grid_tf, H_tf, noise_var, tx_bits, params, mode)
% OFDM receiver with perfect CSI.
%   mode 'onetap' (default): per-subcarrier MMSE, x = conj(h).*y./(|h|^2+s2)
%                            (conventional OFDM; ICI is left as noise)
%   mode 'mmse'            : full-matrix MMSE (ICI-aware, upper bound)
% tx_bits is used ONLY to count errors.

    M = params.M; N = params.N; MN = M * N;
    if nargin < 6, mode = 'onetap'; end

    y = rx_grid_tf(:);
    switch lower(mode)
        case 'onetap'
            d = diag(H_tf);
            x_hat = conj(d) .* y ./ (abs(d).^2 + noise_var);
        case 'mmse'
            x_hat = (H_tf' * H_tf + noise_var * eye(MN)) \ (H_tf' * y);
        otherwise
            error('Unknown mode %s', mode);
    end
    rx_syms = reshape(x_hat, M, N);

    rx_bits = zeros(2*MN, 1);
    rx_bits(1:2:end) = double(real(x_hat) < 0);
    rx_bits(2:2:end) = double(imag(x_hat) < 0);

    num_errors = sum(tx_bits(:) ~= rx_bits(:));
    ber = num_errors / length(tx_bits);
end

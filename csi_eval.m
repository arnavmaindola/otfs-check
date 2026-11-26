function e = csi_eval(y_o, H_o, y_f, H_f, nv, bits, p, s_o, s_f, s_d, E_o, E_f, E_d)
% CSI_EVAL  Bit errors [OTFS-MMSE, OFDM one-tap, OFDM full MMSE] with H_est = H + error.
%   s_* : error amplitude relative to sqrt(mean power): s_o (OTFS matrix),
%   s_f (OFDM matrix), s_d (one-tap diagonal, rel. to mean subcarrier-gain power).
    MN = size(H_o, 1);
    He_o = H_o + s_o*sqrt(mean(abs(H_o(:)).^2))*E_o;
    He_f = H_f + s_f*sqrt(mean(abs(H_f(:)).^2))*E_f;
    He_1 = H_f + diag(s_d*sqrt(mean(abs(diag(H_f)).^2))*E_d);
    e = zeros(1,3);
    [~,~,~,e(1)] = otfs_receiver(y_o, He_o, nv, bits, p);
    [~,~,~,e(2)] = ofdm_receiver(y_f, He_1, nv, bits, p, 'onetap');
    [~,~,~,e(3)] = ofdm_receiver(y_f, He_f, nv, bits, p, 'mmse');
    %#ok<NASGU>
end

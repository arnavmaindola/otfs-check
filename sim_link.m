function ber = sim_link(params, snr_list, num_frames, seed, do_ofdm)
% SIM_LINK  Monte-Carlo BER of OTFS (and optionally OFDM) over the same
% channel realisations.
%   ber(:,1) OTFS + MMSE        ber(:,2) OFDM one-tap MMSE
%   ber(:,3) OFDM full MMSE (ICI-aware)
% Frame f uses rng(seed+f): bits, channel and noise are identical across
% velocities/models (common random numbers => smooth curves, few frames).

    if nargin < 5, do_ofdm = true; end
    ns = numel(snr_list);
    err = zeros(ns, 3);
    nbits = 2 * params.M * params.N;

    for f = 1:num_frames
        rng(seed + f);
        chan = channel_state(params);
        p = params; p.cp_len = chan.max_delay;
        bits = randi([0 1], nbits, 1);

        [tx_otfs, ~, ~] = otfs_transmitter(p, bits);
        if do_ofdm, tx_ofdm = ofdm_transmitter(p, bits); end

        for s = 1:ns
            p.snr_db = snr_list(s);
            [y, H, nv] = otfs_channel(tx_otfs, p, chan);
            [~, ~, ~, e] = otfs_receiver(y, H, nv, bits, p);
            err(s, 1) = err(s, 1) + e;

            if do_ofdm
                [y, H, nv] = ofdm_channel(tx_ofdm, p, chan);
                [~, ~, ~, e1] = ofdm_receiver(y, H, nv, bits, p, 'onetap');
                [~, ~, ~, e2] = ofdm_receiver(y, H, nv, bits, p, 'mmse');
                err(s, 2) = err(s, 2) + e1;
                err(s, 3) = err(s, 3) + e2;
            end
        end
    end
    ber = err / (nbits * num_frames);
end

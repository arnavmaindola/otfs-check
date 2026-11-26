function chan = channel_state(params)
% CHANNEL_STATE  Draws ONE random realisation of the EVA multipath channel.
% Both OTFS and OFDM links are driven from the same 'chan', so the
% comparison is fair (same fades, same Dopplers).
%
% params: M, delta_f, fc, velocity_kmh, fading_model
% chan  : lt (integer delay taps), h (P x 1 complex gains),
%         nu (P x 1 Doppler shifts [Hz]), Ts, max_delay

    M = params.M;
    Ts = 1 / (M * params.delta_f);

    % Standard 3GPP EVA profile (9 taps)
    tau_sec = [0, 30, 150, 310, 370, 710, 1090, 1730, 2510] * 1e-9;
    p_db    = [0.0, -1.5, -1.4, -3.6, -0.6, -9.1, -7.0, -12.0, -16.9];
    p = 10.^(p_db / 10);
    p = (p / sum(p)).';            % unit total average power
    P = numel(p);

    lt = round(tau_sec / Ts).';    % integer delay taps (sample spacing Ts)

    % Jakes Doppler per path: nu = nu_max*cos(theta)
    nu_max = (params.velocity_kmh / 3.6) * params.fc / 3e8;
    use_crn = isfield(params, 'crn');   % optional: externally supplied variates (see crn_variates)
    if use_crn, nu = nu_max * cos(2*pi*params.crn.ang);
    else,       nu = nu_max * cos(2*pi*rand(P, 1)); end

    cn = @() randn(P, 1) + 1i*randn(P, 1);
    ph = @() exp(1i*2*pi*rand(P, 1));

    % Every model is normalised so that E|h_p|^2 = p_p
    if use_crn
        h = crn_gains(params, p, P);
    else
    switch lower(params.fading_model)
        case 'rayleigh'
            h = sqrt(p/2) .* cn();

        case 'rician'      % K = 6 dB of total power in one LOS ray on tap 1
            K = 10^(6/10);
            h = sqrt(p/(K+1)/2) .* cn();
            h(1) = h(1) + sqrt(K/(K+1)) * exp(1i*2*pi*rand());

        case 'rician-all-taps'   % K = 6 dB (params.K_db) on EVERY tap; = kappa-mu with mu = 1
            K = 10^(6/10);
            if isfield(params,'K_db'), K = 10^(params.K_db/10); end
            h = sqrt(p/(K+1)/2) .* cn() + sqrt(p*K/(K+1)) .* ph();

        case 'nakagami-m'  % m = 2
            m = 2;
            h = sqrt(p/m .* gamma_rnd(m, P)) .* ph();

        case 'alpha-mu'    % alpha = 1.2, mu = 0.8 (severe fading)
            alpha = 1.2; mu = 0.8;
            if isfield(params,'alpha'), alpha = params.alpha; end   % optional overrides
            if isfield(params,'mu'),    mu    = params.mu;    end
            r_hat = sqrt(p * mu^(2/alpha) * gamma(mu) / gamma(mu + 2/alpha));
            h = r_hat .* (gamma_rnd(mu, P) / mu).^(1/alpha) .* ph();

        case 'kappa-mu'    % kappa = 3, mu = 2
            kappa = 3; mu = 2;
            if isfield(params,'kappa'), kappa = params.kappa; end   % optional overrides
            if isfield(params,'mu'),    mu    = params.mu;    end   % integer mu
            d  = kappa * p / (1 + kappa);
            s2 = p / (2*mu*(1 + kappa));
            Gc = sqrt(d / (2*mu)).' + sqrt(s2).' .* randn(2*mu, P);
            h  = sqrt(sum(Gc.^2, 1)).' .* ph();

        otherwise
            error('Unknown fading model: %s', params.fading_model);
    end
    end

    chan.lt = lt; chan.h = h; chan.nu = nu; chan.Ts = Ts;
    chan.max_delay = max(lt);
end

function h = crn_gains(params, p, P)
% Gains for every model built from the SAME externally drawn variates:
%   crn.n (8 x P) N(0,1), crn.ph (P x 1) U(0,1), crn.u (P x 1) U(0,1).
% Rayleigh/Rician use n(1:2,:); kappa-mu uses n(1:2*mu,:); alpha-mu uses u.
    V = params.crn; ph = exp(1i*2*pi*V.ph);
    cn = (V.n(1,:) + 1i*V.n(2,:)).';
    switch lower(params.fading_model)
        case 'rayleigh'
            h = sqrt(p/2) .* cn;
        case 'rician'
            K = 10^(6/10); if isfield(params,'K_db'), K = 10^(params.K_db/10); end
            h = sqrt(p/(K+1)/2) .* cn;
            h(1) = h(1) + sqrt(K/(K+1)) * ph(1);
        case 'rician-all-taps'
            K = 10^(6/10); if isfield(params,'K_db'), K = 10^(params.K_db/10); end
            h = sqrt(p/(K+1)/2) .* cn + sqrt(p*K/(K+1)) .* ph;
        case 'nakagami-m'      % m = 2: Gamma(2) = sum of 4 squared N(0,1) / 2
            m = 2;
            h = sqrt(p/m .* (sum(V.n(1:4,:).^2, 1).' / 2)) .* ph;
        case 'kappa-mu'
            kappa = params.kappa; mu = params.mu;
            d = kappa*p/(1+kappa); s2 = p/(2*mu*(1+kappa));
            Gc = sqrt(d/(2*mu)).' + sqrt(s2).' .* V.n(1:2*mu, :);
            h = sqrt(sum(Gc.^2, 1)).' .* ph;
        case 'alpha-mu'
            alpha = params.alpha; mu = params.mu;
            r_hat = sqrt(p * mu^(2/alpha) * gamma(mu) / gamma(mu + 2/alpha));
            h = r_hat .* (gammaincinv(V.u, mu) / mu).^(1/alpha) .* ph;
        otherwise
            error('Model %s not supported with crn', params.fading_model);
    end
end

function g = gamma_rnd(k, n)
% Gamma(shape k, scale 1) samples, Marsaglia-Tsang (no toolbox needed).
    g = zeros(n, 1);
    for i = 1:n
        kk = k + (k < 1);
        d = kk - 1/3; c = 1/sqrt(9*d);
        while true
            x = randn(); v = (1 + c*x)^3;
            if v > 0 && log(rand()) < 0.5*x^2 + d - d*v + d*log(v)
                break;
            end
        end
        g(i) = d * v;
        if k < 1
            g(i) = g(i) * rand()^(1/k);
        end
    end
end

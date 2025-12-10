function [s] = generateNLFM(simParams)

% Number of samples per pulse
nPulse = length(simParams.t);

% Amplitude window
U = ones(nPulse, 1);

% Spectral Window
G = abs(simParams.spectral_template);

% Form nonlinear frequency modulated pulse
s = runNLFM(U, G);

% Normalize pulse to unit energy
s = s./sqrt(sum(abs(s).^2));

function [s] = runNLFM(U, G)

    % Principle of stationary phase implementation (based on Fowle 1964)

    % Make column vectors, normalize
    U = U(:)/mean(U);
    G = G(:)/mean(G);

    % Compute Time, Frequency Integrals
    P_t = cumtrapz(U);          % Time integral (left side of (18))
    P_t = P_t/P_t(end);         % Normalize to unit energy
    Q_f = cumtrapz(G);          % Frequency integral (left side of (18))
    Q_f = Q_f/Q_f(end);         % Normalize to unit energy
    
    % Make frequency vector; f = [-0.5 0.5]
    nfft = length(G);
    f = fftshift((0:nfft-1).'/nfft); f(f>=0.5) = f(f>=0.5)-1;

    % Edge Case (Redundant values for notched power spectrum)
    [Q2_f, ii] = unique(Q_f, 'last');
    jj = find(Q_f == Q2_f(end), 1 );
    ii(end) = jj;
    f2 = f(ii);
    
    % Determine Instantaneous Frequency F
    F = interp1(Q2_f,f2,P_t);   % Equation 27 Q^{-1}{P(t)}
    
    % Calculate baseband signal (Equation 28)
    s = U.*exp(1j*2*pi*cumtrapz(F));

 
return
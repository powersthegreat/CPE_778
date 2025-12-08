function [s] = generateCELSI(simParams, s_init)

% Number of samples per pulse
nPulse = length(simParams.t);

% Form binary spectral mask
G = abs(simParams.spectral_template);
wf = fftshift(G <= 0.4);

% figure("Name", "Debuggign PSD Template");
% plot(db(G + 0.001)); hold on;
% plot(db(wf + 0.001));

% Percent power scaling for inequality constraint
gamma = (sum(G<=0)/length(G))/(nPulse)*(10^(6/10));

% Parse initialization waveform
if nargin < 2 || isempty(s_init)
    % No initial waveform provided -> default to LFM
    s_init = generateLFM(simParams);
end
s_init = [real(s_init); imag(s_init)];

% Define autocorrelation mask
wsl = computeCorrelationMask(s_init, simParams.nfft);

% Define optimization parameters
optParams = struct();
optParams.iter = 50;
optParams.p = 8;
optParams.gamma = gamma;
optParams.nPulse = nPulse;
optParams.nfft = simParams.nfft;
optParams.wsl = wsl;
optParams.wf = wf;

% Objective and constraint function definitions
J = @(x) Jx(x, optParams);
C = @(x) Cx(x, optParams);

% % debugging objective and objective gradient
% J = @(x) Jx(x, optParams);
% VJ = @(x) VJx(x, optParams);
% checkGradient(J, VJ, s_init);

% % debugging inequality constraint and inequality constraint gradient
% G = @(x) Gx(x, optParams);
% VG = @(x) VGx(x, optParams);
% checkGradient(G, VG, s_init);

% % debugging equality constraint and equality constraint gradient
% H = @(x) Hx(x, optParams);
% VH = @(x) VHx(x, optParams);
% checkGradient(H, VH, s_init);

% Run optimization usign fmincon toolbox
options = optimoptions('fmincon','MaxFunctionEvaluations',1e10,'MaxIterations',optParams.iter, ...
    'OptimalityTolerance',1e-12, 'StepTolerance',1e-12,'Algorithm','interior-point',...
    'SpecifyObjectiveGradient', false, 'SpecifyConstraintGradient', false, 'Display', 'iter');
x_opt = fmincon(J, s_init, [],[],[],[],[],[], C, options);

% Reconstruct waveform from real and imaginary
s = x_opt(1:nPulse) + 1i*x_opt(nPulse+1:end);

% Normalize pulse to unit energy
s = s./sqrt(sum(abs(s).^2));



%% Objective and Constraint Functions

function [J, VJ] = Jx(x, optParams)
    
    nPulse = optParams.nPulse;
    nfft = optParams.nfft;
    wsl = optParams.wsl;
    p = optParams.p;

    % Reconstruct waveform from real and imaginary
    s = x(1:nPulse) + 1i*x(nPulse+1:end);

    % Evaluate cost function
    ccorr = ifft( fft(s, nfft) .* conj(fft(s, nfft)));
    J = log(sum(abs(wsl.*ccorr).^p));

    % Call objective gradient function
    VJ = VJx(x, optParams);

    % % Debugging ACF
    % figure("Name", "Debugging ACF Mask");
    % plot(fftshift(db(ccorr))); hold on;
    % plot(fftshift(db(wsl + 0.001)));
    % ylim([-75 5]);
    
return

function [VJ] = VJx(x, optParams)

    nfft = optParams.nfft;
    wsl = optParams.wsl;
    p = optParams.p;
    nPulse = optParams.nPulse;

    % Reconstruct waveform from real and imaginary
    s = x(1:nPulse) + 1i.*x(nPulse+1:end);

    % Evaluate cost function
    ccorr = ifft(fft(s, nfft) .* conj(fft(s, nfft)));
    Jsl = sum(abs(wsl.*ccorr).^p);

    % Redundant pre-calculations
    temp = (abs(ccorr).^(p-2)).*ccorr;
    tempsl = fft(wsl.*temp);

    % CELSI gradient
    S = fft(s, nfft);
    c = ifft(S .* tempsl);
    g_s = p * c(1:nPulse);
    g_re = 2 * real(g_s) / Jsl;
    g_im = 2 * imag(g_s) / Jsl;
    VJ = [g_re; g_im];

return

function [G] = Gx(x, optParams)

    wf = optParams.wf;
    gamma = optParams.gamma;
    nfft = optParams.nfft;
    nPulse = optParams.nPulse;

    % Reconstruct waveform from real and imaginary
    s = x(1:nPulse) + 1i.*x(nPulse+1:end);

    % Inequality constraint function
    X = fft(s, nfft) / sqrt(nfft);
    G = sum(wf .* abs(X).^2) - gamma * sum(abs(s).^2);
    
return

function [VG] = VGx(x, optParams)

    wf     = optParams.wf;
    gamma  = optParams.gamma;
    nfft   = optParams.nfft;
    nPulse = optParams.nPulse;

    % Reconstruct waveform
    s = x(1:nPulse) + 1i*x(nPulse+1:end);

    % FFT with 1/sqrt(nfft)
    X = fft(s, nfft) / sqrt(nfft);

    % === CORRECT GRADIENT: include 1/sqrt(nfft) ===
    g_s_full = ifft(wf .* X) / sqrt(nfft);

    % Truncate to N time-domain samples and subtract energy term
    g_s = g_s_full(1:nPulse) - gamma * s;

    % Convert to real + imag variables
    VG_real = 2 * real(g_s);
    VG_imag = 2 * imag(g_s);
    VG = [VG_real; VG_imag];
return

function H = Hx(x,optParams)
    n = optParams.nPulse;
    s = x(1:n) + 1i*x(n+1:end);
    H = abs(s).^2 - 1;
return

function VH = VHx(x,optParams)

    n = optParams.nPulse;
    s = x(1:n) + 1i*x(n+1:end);

    VH = zeros(n,2*n);
    VH(:,1:n)     = diag(2*real(s));
    VH(:,n+1:end) = diag(2*imag(s));
    VH = VH.';

return

function [G, H, VG, VH] = Cx(x, optParams)

    G = Gx(x, optParams);
    VG = VGx(x, optParams);
    H = Hx(x, optParams);
    VH = VHx(x, optParams);

return


function [s] = generateCELSI_base(simParams, s_init)

% Number of samples per pulse
nPulse = length(simParams.t);

% Form binary spectral mask
G = abs(simParams.spectral_template);
wf = fftshift(G <= 0);

% Percent power scaling for inequality constraint
gamma = (sum(G<=0)/length(G))/(nPulse)*(10^(6/10));

% Parse initialization waveform
if nargin < 2 || isempty(s_init)
    % No initial waveform provided -> default to LFM
    s_init = generateLFM(simParams);
end
phi_init = angle(s_init);

% Define autocorrelation mask
wsl = computeCorrelationMask(s_init, simParams.nfft);

% Define optimization parameters
optParams = struct();
optParams.iter = 100;
optParams.p = 8;
optParams.gamma = gamma;
optParams.nPulse = nPulse;
optParams.nfft = simParams.nfft;
optParams.wsl = wsl;
optParams.wf = wf;

% % Objective and constraint function definitions
% J = @(x) Jx(x, optParams);
% C = @(x) Cx(x, optParams);

% debugging gradient with finite difference
J = @(x) Jx(x, optParams);
VJ = @(x) VJx(x, optParams);

checkGradient(J, VJ, s_init);
s = 0;

% % Run using optimization toolbox
% options = optimoptions('fmincon','MaxFunctionEvaluations',1e10,'MaxIterations',optParams.iter, ...
%     'OptimalityTolerance',1e-10, 'StepTolerance',1e-10,'Algorithm','interior-point',...
%     'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true, ...
%     'Display', 'iter');
% phi_opt = fmincon(J, phi_init, [],[],[],[],[],[], C, options);

% % Form signal from optimized phase history
% s = exp(1i*phi_opt);

% % Normalize pulse to unit energy
% s = s./sqrt(sum(abs(s).^2));


%% Objective and Constraint Functions

function [J, VJ] = Jx(x, optParams)
    
    nPulse = optParams.nPulse;
    nfft = optParams.nfft;
    wsl = optParams.wsl;
    p = optParams.p;

    % Create signal
    s = exp(1i*x)./sqrt(nPulse);

    % Evaluate cost function
    ccorr = ifft( fft(s, nfft) .* conj(fft(s, nfft)));
    J = log(sum(abs(wsl.*ccorr).^p));
    VJ = VJx(x, optParams);

    % % Debugging ACF
    % figure("Name", "Debugging ACF Mask");
    % plot(db(ccorr)); hold on;
    % plot(db(wsl));
    
return

function [VJ] = VJx(x, optParams)

    nPulse = optParams.nPulse;
    nfft = optParams.nfft;
    wsl = optParams.wsl;
    p = optParams.p;
    
    % Create signal
    s = exp(1i*x)./sqrt(nPulse);

    % Evaluate cost function
    ccorr = ifft( fft(s, nfft) .* conj(fft(s, nfft)));
    Jsl = sum(abs(wsl.*ccorr).^p);

    % Redundant pre-calculations
    temp = (abs(ccorr).^(p-2)).*ccorr;
    tempsl = fft(wsl.*temp);

    % CELSI gradient
    c = ifft(fft(s, nfft) .* tempsl);
    g = p*imag(conj(s).*c(1:size(s,1),:));
    VJ = 2*g/Jsl;

return

function [G] = Gx(x, optParams)

    wf = optParams.wf;
    gamma = optParams.gamma;
    nfft = optParams.nfft;
    nPulse = optParams.nPulse;
    
    s = exp(1i*x)./sqrt(nPulse);
    G = sum(wf.*((abs(fft(s,nfft)/sqrt(nfft)).^2))) - gamma*sum(abs(s).^2);

    % figure("Name", "Debugging PSD Mask");
    % plot(db(abs(fft(s,nfft)/sqrt(nfft)).^2)); hold on;
    % plot(db(wf + 0.001));
    
return

function [VG] = VGx(x, optParams)

    wf = optParams.wf;
    gamma = optParams.gamma;
    nfft = optParams.nfft;
    nPulse = optParams.nPulse;
    
    s = exp(1i*x)./sqrt(nPulse);
    temp = ifft(fft(s, nfft).*wf);
    VG   = 2*imag(conj(s).*(temp(1:nPulse)-gamma*s));
   
return

function [G, H, VG, VH] = Cx(x, optParams)

    G = Gx(x, optParams);
    VG = VGx(x, optParams);
    H = []; VH = [];

return


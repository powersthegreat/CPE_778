function [s] = generateLFM(simParams)

    % Start and stop frequencies
    B1 = -simParams.b/2; B2 = simParams.b/2;

    % Local time support
    nPulse = length(simParams.t);
    t = 0:nPulse-1;

    % Generate linear frequency modulated pulse
    s = exp(2i*pi*(B1*t+.5*(B2-B1)/nPulse*t.^2)).';

    % Normalize pulse to unit energy
    s = s./sqrt(sum(abs(s).^2));

end

function [G] = Gx(x, optParams)

    wf    = optParams.wf;
    gamma = optParams.gamma;
    nfft  = optParams.nfft;
    nPulse = optParams.nPulse;

    % Reconstruct waveform from real and imaginary parts
    s = x(1:nPulse) + 1i*x(nPulse+1:end);

    % FFT
    X = fft(s, nfft) / sqrt(nfft);

    % Inequality constraint
    G = sum(wf .* abs(X).^2) - gamma * sum(abs(s).^2);
end

function [VG] = VGx(x, optParams)

    wf    = optParams.wf;
    gamma = optParams.gamma;
    nfft  = optParams.nfft;
    nPulse = optParams.nPulse;

    % Reconstruct waveform
    s = x(1:nPulse) + 1i*x(nPulse+1:end);

    % FFT(s)
    X = fft(s, nfft) / sqrt(nfft);

    % === Gradient in the Fourier domain ===
    % g_s = IFFT(w .* X) - gamma * s
    g_s = ifft(wf .* X) - gamma * s;

    % === Convert to gradient wrt real + imag variables ===
    VG_real = 2 * real(g_s);
    VG_imag = 2 * imag(g_s);

    VG = [VG_real; VG_imag];
end

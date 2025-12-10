function [mf_out, psd_in, bw_6dB] = plotWaveform(waveform, sample_rate)
    %   Generates a single inline (1x3) figure showing:
    %     (1) Correlation Response
    %     (2) Power Spectral Density & 6-dB Bandwidth
    %     (3) Temporal Magnitude (Envelope)


    waveform1 = waveform(:);
    waveform2 = waveform1;

    N = length(waveform1);
    nfft = 2^nextpow2(2*N);

    time_axis = (-nfft/2:nfft/2-1) / sample_rate * 1e6;   % µs
    freq_axis = (-nfft/2:nfft/2-1) * (sample_rate/nfft); % Hz

    % (1) Matched Filter / Correlation Response
    mf_out = abs(fftshift(ifft(fft(waveform1,nfft) .* conj(fft(waveform2,nfft)))));
    mf_norm = 20*log10(mf_out / max(mf_out));

    % (2) Power Spectrum + 6-dB Bandwidth
    psd_in = fftshift(abs(fft(waveform1, nfft)).^2);
    psd_norm = 10*log10(psd_in / max(psd_in));

    % indices above -6 dB
    idx_above = find(psd_norm >= -6);
    if isempty(idx_above)
        f_low = NaN; f_high = NaN; bw_6dB = NaN;
    else
        f_low  = freq_axis(idx_above(1));
        f_high = freq_axis(idx_above(end));
        bw_6dB = f_high - f_low;
    end

    % (3) Magnitude Envelope
    waveform_mag = abs(waveform1);
    time_mag = (0:N-1) / sample_rate * 1e6;

    % Plotting
    figure('Name','Waveform Analysis (1x3)', 'Color','w');
    tiledlayout(1,3, 'TileSpacing','compact', 'Padding','compact');

    nexttile;
    plot(time_axis, mf_norm, 'LineWidth', 1.5);
    xlabel('Time (µs)');
    ylabel('Correlation (dB)');
    title('Correlation Response');
    ylim([-75 3]);
    grid on;

    nexttile;
    plot(freq_axis/1e6, psd_norm, 'LineWidth', 1.5); hold on;
    xlabel('Frequency (MHz)');
    ylabel('PSD (dB)');
    title(sprintf('Power Spectrum'));
    xlim([-225 225]);
    ylim([-75 3]);
    grid on;

    nexttile;
    plot(time_mag, waveform_mag, 'LineWidth', 1.5);
    xlabel('Time (µs)');
    ylabel('Amplitude');
    title('Waveform Magnitude');
    ylim([min(waveform_mag)-0.05, max(waveform_mag)+0.05]);
    grid on;

end

function [psd_in, bw_6dB] = plotPowerSpectrum(waveform, sample_rate)
% plotPowerSpectrum  Plot the normalized power spectral density of a waveform
% and compute its 6-dB bandwidth.
%
%   [psd_in, bw_6dB] = plotPowerSpectrum(waveform, sample_rate)
%
%   Inputs:
%       waveform     : Input signal vector
%       sample_rate  : Sampling rate in Hz
%
%   Outputs:
%       psd_in       : Power spectral density (unnormalized)
%       bw_6dB       : 6-dB bandwidth in Hz
%
%   This function computes the FFT-based power spectral density (PSD),
%   plots the normalized response in dB versus frequency, and marks
%   the -6 dB bandwidth limits.

    % Ensure column vector
    waveform = waveform(:);

    % Zero-pad for better frequency resolution
    nfft = 2^nextpow2(2 * length(waveform));

    % Frequency axis (centered around zero)
    freq_axis = (-nfft/2:nfft/2-1) * (sample_rate/nfft);

    % Power spectrum (normalized)
    psd_in = fftshift(abs(fft(waveform, nfft)).^2);
    psd_norm = 10*log10(psd_in / max(psd_in));

    % === Compute 6-dB Bandwidth ===
    % Find where the PSD crosses -6 dB
    idx_above = find(psd_norm >= -6);
    if isempty(idx_above)
        warning('No points found above -6 dB threshold.');
        bw_6dB = NaN;
    else
        f_low  = freq_axis(idx_above(1));
        f_high = freq_axis(idx_above(end));
        bw_6dB = f_high - f_low;
    end

    % === Plot ===
    figure('Name', 'Power Spectral Response');
    plot(freq_axis/1e6, psd_norm, 'LineWidth', 1.5); hold on;
    xlabel('Frequency (MHz)');
    ylabel('Normalized Power (dB)');
    title('Waveform Power Spectral Response');
    ylim([-75 0]);
    % xlim([-200 200]);
    grid on;

    % Mark -6 dB points
    ylims = ylim;
    plot([f_low f_low]/1e6, ylims, 'r--', 'LineWidth', 1.2);
    plot([f_high f_high]/1e6, ylims, 'r--', 'LineWidth', 1.2);
    text(f_low/1e6, ylims(1)+3, sprintf('%.3f MHz', f_low/1e6), 'Color', 'r', 'HorizontalAlignment', 'right');
    text(f_high/1e6, ylims(1)+3, sprintf('%.3f MHz', f_high/1e6), 'Color', 'r', 'HorizontalAlignment', 'left');

    % Annotate bandwidth
    bw_MHz = bw_6dB / 1e6;
    legend('PSD (dB)', sprintf('6-dB BW = %.3f MHz', bw_MHz), 'Location', 'best');
    hold off;

    % Display numeric output
    fprintf('6-dB Bandwidth: %.3f MHz (%.0f Hz)\n', bw_MHz, bw_6dB);

end

function [mf_out] = plotCorrelationResponse(waveform1, waveform2, sample_rate)
% plotCorrelationResponse  Compute and plot the correlation response
%   between two equally sized waveforms
%
%   plotCorrelationResponse(waveform1, waveform2, sample_rate)
%
%   Inputs:
%       waveform1    : First input signal vector
%       waveform2    : Second input signal vector (same length as waveform1)
%       sample_rate  : Sampling rate in Hz
%
%   This function computes the matched filter/correlation response
%   using FFT-based convolution and plots it in dB versus time.

    % Ensure column vectors
    waveform1 = waveform1(:);
    waveform2 = waveform2(:);

    if length(waveform1) ~= length(waveform2)
        error('Waveform inputs must be the same length.');
    end

    % Zero-pad for correlation
    nfft = 2^nextpow2(2 * length(waveform1));

    % Matched filter (correlation)
    mf_out = abs(fftshift(ifft(fft(waveform1, nfft) .* conj(fft(waveform2, nfft)))));

    % Normalize and convert to dB
    mf_norm = 20*log10(mf_out / max(mf_out));

    % Time axis
    time_axis = (-nfft/2:nfft/2-1) / sample_rate;

    % Plot
    figure('Name', 'Correlation Response');
    plot(time_axis*1e6, mf_norm, 'LineWidth', 1.5);
    xlabel('Time (µs)');
    ylabel('Normalized Correlation (dB)');
    title('Waveform Correlation Response');
    ylim([-75 0]);
    grid on;
end

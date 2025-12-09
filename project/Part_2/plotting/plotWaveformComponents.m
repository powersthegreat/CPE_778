function plotWaveformComponents(waveform, sample_rate)
% plotWaveformComponents  Plot real, imaginary, magnitude, and instantaneous frequency
%
%   plotWaveformComponents(waveform, sample_rate)
%
%   Inputs:
%       waveform     : Input complex signal vector
%       sample_rate  : Sampling rate in Hz
%
%   This function plots four subplots stacked vertically:
%       - Real part of the waveform
%       - Imaginary part of the waveform
%       - Magnitude (absolute value) of the waveform
%       - Instantaneous frequency (Hz)

    % Ensure column vector
    waveform = waveform(:);

    % Time axis in microseconds
    time_axis = (0:length(waveform)-1) / sample_rate * 1e6;

    % Create figure
    figure('Name', 'Waveform Components');

    % Magnitude
    power = sum(abs(waveform).^2);
    plot(time_axis, abs(waveform), 'LineWidth', 1.5);
    xlabel('Time (µs)');
    ylabel('Amplitude');
    title(sprintf('Magnitude (Normalized Power = %.2f)', power));
    ylim([-.1, 1.1]);
    grid on;

end

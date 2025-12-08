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

    % Compute instantaneous frequency
    phase_unwrapped = unwrap(angle(waveform));
    inst_freq = diff(phase_unwrapped) * sample_rate / (2*pi);
    inst_phase = 2 * pi * cumtrapz(time_axis(2:end), inst_freq);

    % Create figure
    figure('Name', 'Waveform Components', 'Position', [100 100 1200 800]);

    % Real part
    subplot(2, 2, 1);
    plot(time_axis, real(waveform), 'b', 'LineWidth', 1.5); hold on;
    plot(time_axis, imag(waveform), 'r', 'LineWidth', 1.5);
    xlabel('Time (µs)');
    ylabel('Amplitude');
    title('Real/Imag Components');
    grid on;

    % Magnitude
    power = sum(abs(waveform).^2);
    subplot(2, 2, 2);
    plot(time_axis, abs(waveform), 'LineWidth', 1.5);
    xlabel('Time (µs)');
    ylabel('Amplitude');
    title(sprintf('Magnitude (Power = %.2f)', power));
    ylim([-.1, 1.1]);
    grid on;

    % Plot instantaneous frequency
    subplot(2, 2, 3);
    plot(time_axis(2:end), inst_freq, 'LineWidth', 1.5);
    xlabel('Time (µs)');
    ylabel('Frequency (Hz)');
    title('Instantaneous Frequency');
    grid on;

    % Plot instantaneous phase progression
    subplot(2, 2, 4);
    plot(time_axis(2:end), unwrap(inst_phase), 'b', 'LineWidth', 1.5);
    xlabel('Time (µs)');
    ylabel('Phase (rad)');
    title('Instantaneous Phase Progression');
    grid on;


end

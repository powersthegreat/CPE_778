% clc; clear;
addpath('plotting');
addpath('helpers');



%% Define Design and Simulation Parameters

% define design parameters
simParams = struct();
simParams.bandwidth       = 100e6;
simParams.sample_rate     = .5e9;
simParams.pulse_width     = 1e-6;

% set simulation paramters
simParams.b               = simParams.bandwidth / simParams.sample_rate;
simParams.t               = ((1/simParams.sample_rate):(1/simParams.sample_rate):simParams.pulse_width).';
simParams.nfft            = 2*(simParams.pulse_width*simParams.sample_rate)-1;

% generate spectral template
simParams.spectral_template = windowKaiser(simParams.nfft, 0.585*(1/simParams.b), -80);



%% Generate Waveform

% generate LFM waveform
waveform = generateLFM(simParams);

% % generate NLFM waveform
% waveform = generateNLFM(simParams);

plotWaveformComponents(waveform / max(waveform), simParams.sample_rate);
plotPowerSpectrum(waveform, simParams.sample_rate);
plotCorrelationResponse(waveform, waveform, simParams.sample_rate);

% generate CELSI waveform
waveform = generateCELSI(simParams);



%% Plot Outputs

plotWaveformComponents(waveform / max(waveform), simParams.sample_rate);
plotPowerSpectrum(waveform, simParams.sample_rate);
plotCorrelationResponse(waveform, waveform, simParams.sample_rate);
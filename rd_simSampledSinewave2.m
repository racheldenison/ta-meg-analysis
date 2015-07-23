% rd_simSampledSinewave2.m

%% setup
sinewave = @(A,f,t,ph) A*sin(2*pi*f*t + ph);

dur = 2;

refRate = 60;
tStim = -0.5:1/refRate:dur;
stimAmp = 1;
stimFreq = 15;
stimFreq2 = 20;
phase = pi/2;
phase2 = 0;

Fs = 1000;
tResponse = -0.5:1/Fs:dur;
noiseStd = 0;

plotFigs = 1;

%% STIMULUS
%% constructed neural response - use 120 Hz
% assumes that neurons respond to both onsets and offsets
% so if a stimulus at 60 Hz is [1 0 0], the neural response would be 
% [1 1 0]. sampled at 120 Hz, this would be [1 0 1 0 0 0].
% % 30 Hz
% stim = repmat([1 0 0 0],1,301/4);
% % 40 Hz
% stim = repmat([1 0 1 0 0 0],1,301/6);
% stim(301) = 1;

%% make stimulus time series
stim = sinewave(stimAmp,stimFreq,tStim,phase);
if ~isempty(stimFreq2)
    stim = stim + sinewave(stimAmp,stimFreq2,tStim,phase2);
end
stim(tStim<0) = 0; % zero baseline

if plotFigs % && 0
    f(1) = figure;
    plot(tStim, stim, '.-')
    xlabel('time (s)')
    ylabel('amp')
    title(sprintf('%.2f Hz', stimFreq))
end

%% wavelet on stim
width = 10;
[spectrum,freqoi,timeoi] = ft_specest_wavelet(stim, tStim, 'freqoi', stimFreq, 'width', width);
specAmp = abs(squeeze(spectrum))';

if plotFigs
    figure
    plot(tStim, specAmp)
    xlabel('time (ms)')
    ylabel('spec amp')
    title('wavelet')
end

%% time-frequency on stim
taper          = 'hanning';
foi            = 1:50;
t_ftimwin      = 10 ./ foi;
toi            = tStim;
tfAmps = [];
[spectrum,ntaper,freqoi,timeoi] = ft_specest_mtmconvol(stim, tStim, ...
    'timeoi', toi, 'freqoi', foi, 'timwin', t_ftimwin, ...
    'taper', taper, 'dimord', 'chan_time_freqtap');
specAmp = abs(squeeze(spectrum))';

if plotFigs
    ytick = 10:10:numel(foi);
    xtick = 1:30:numel(toi);
    figure
    imagesc(specAmp)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick);
    title('mtmconvol')
end

%% SIMULATED NEURAL RESPONSE
%% make neural response 
% by assuming that neurons respond to differences
stimUp = resample(stim, Fs, refRate);
ssvef = abs(diff(stimUp));
% ssvef = resample(ssvef, Fs, refRate);
if length(ssvef) > length(tResponse)
    ssvef = ssvef(1:length(tResponse));
elseif length(tResponse) > length(ssvef)
    tResponse = tResponse(1:length(ssvef));
end
noise = noiseStd.*randn(size(tResponse));

response = ssvef + noise;

if plotFigs
    scaleFactor = max(stim)./max(ssvef);
    figure
    hold on
    plot(tStim, stim, '.-')
    plot(tResponse, ssvef*scaleFactor, 'r')
    xlabel('time (s)')
    ylabel('amp')
    title(sprintf('%.2f Hz', stimFreq))
end

%% wavelet on simulated neural response
[spectrum,freqoi,timeoi] = ft_specest_wavelet(response, tResponse, 'freqoi', stimFreq);
specAmp = abs(squeeze(spectrum))';

if plotFigs
    figure
    plot(tResponse, specAmp)
    xlabel('time (ms)')
    ylabel('spec amp')
end

%% time-frequency on simulated neural response
taper          = 'hanning';
foi            = 1:80;
t_ftimwin      = 10 ./ foi;
toi            = tResponse;
tfAmps = [];
[spectrum,ntaper,freqoi,timeoi] = ft_specest_mtmconvol(response, tResponse, ...
    'timeoi', toi, 'freqoi', foi, 'timwin', t_ftimwin, ...
    'taper', taper, 'dimord', 'chan_time_freqtap');
specAmp = abs(squeeze(spectrum))';

if plotFigs
    ytick = 10:10:numel(foi);
    xtick = 1:500:numel(toi);
    figure
    imagesc(specAmp)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick);
end


% rd_simSampledSinewave.m

%% setup
sinewave = @(A,f,t,ph) A*sin(2*pi*f*t + ph);

dur = 2;

refRate = 120;
tStim = -0.5:1/refRate:dur;
stimAmp = 1;
stimFreq = 30;
stimFreq2 = 40;

Fs = 1000;
tResponse = -0.5:1/Fs:dur;
noiseStd = 0;

plotFigs = 1;

%% STIMULUS
%% make stimulus time series
stim = sinewave(stimAmp,stimFreq,tStim,0);
if ~isempty(stimFreq2)
    stim = stim + sinewave(stimAmp,stimFreq2,tStim,0);
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

% foi = 1:50;
% [spectrum,freqoi,timeoi] = ft_specest_wavelet(stim, tStim, 'freqoi', foi);
% specAmp = abs(squeeze(spectrum));
% 
% if plotFigs
%     ytick = 10:10:numel(freqoi);
%     xtick = 51:50:numel(timeoi);
%     figure
%     imagesc(specAmp)
%     rd_timeFreqPlotLabels(timeoi,freqoi,xtick,ytick);
% end

%% mtmconvol on stim
taper          = 'hanning'; % 'dpss'
toi            = tStim;
foi            = stimFreq;
t_ftimwin      = 12 ./ foi;
tapsmofrq      = foi * 0.8; % for dpss 
tfAmps = [];
[spectrum,ntaper,freqoi,timeoi] = ft_specest_mtmconvol(stim, tStim, ...
    'timeoi', toi, 'freqoi', foi, 'timwin', t_ftimwin, ...
    'taper', taper, 'tapsmofrq', tapsmofrq, 'dimord', 'chan_time_freqtap');
specAmp = abs(squeeze(spectrum))';

if plotFigs
    figure
    plot(tStim, specAmp)
    xlabel('time (ms)')
    ylabel('spec amp')
    title('mtmconvol')
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
% by resampling the stimulus -- right thing to do?
ssvef = resample(stim, Fs, refRate);
ssvef = ssvef(1:length(tResponse));
noise = noiseStd.*randn(size(tResponse));

response = ssvef + noise;

if plotFigs
    figure(f(1));
    hold on
    plot(tResponse, ssvef, 'r')
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
foi            = 1:50;
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


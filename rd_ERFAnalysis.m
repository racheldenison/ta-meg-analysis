% rd_ERFAnalysis.m

%% Setup
% exptDir = '/Local/Users/denison/Data/TAPilot/MEG';
exptDir = '/Volumes/DRIVE1/DATA/rachel/MEG/TAPilot/MEG';
sessionDir = 'R0817_20140820';
fileBase = 'R0817_TAPilot_8.20.14';
analStr = 'ebi';
excludeTrialsFt = 1;

dataDir = sprintf('%s/%s', exptDir, sessionDir);

switch analStr
    case ''
        filename = sprintf('%s/%s.sqd', dataDir, fileBase);
        savename = sprintf('%s/mat/%s_erf_workspace.mat', dataDir, fileBase);
        figDir = sprintf('%s/figures/raw', dataDir);
    otherwise
        filename = sprintf('%s/%s_%s.sqd', dataDir, fileBase, analStr);
        savename = sprintf('%s/mat/%s_%s_erf_workspace.mat', dataDir, fileBase, analStr);
        figDir = sprintf('%s/figures/%s', dataDir, analStr);
end
if ~exist(figDir,'dir')
    mkdir(figDir)
end

% trigChan = 160:167;
trigChan = 164:165; % targets
megChannels = 0:156;
channelSets = {0:39,40:79,80:119,120:156};

% channels
switch sessionDir
    case 'R0890_20140806'
        % badChannels = [10 11 115]; % R0890, 48-->49, 150-->152
        badChannels = [10 11 115 49 152]; % R0890
        highSNRChannelsL = [26 60 14 92]; % R0890
        highSNRChannelsR = [1 50 7 8]; % R0890
    case 'R0817_20140820'
        badChannels = [115 152]; % R0817
        highSNRChannelsL = [92 60 15 14]; % R0817
        highSNRChannelsR = [51 1 50 39]; % R0817
    otherwise
        error('sessionDir not found')
end
% badChannels = [];

pickedChannels = {highSNRChannelsL, highSNRChannelsR};
pickedChannelsLabels = {'L','R'};

tstart = -1000; % for targets
tstop = 2000;
t = tstart:tstop;

% trigNames = {'fastL-attL','fastL-attR','fastR-attL','fastR-attR',...
%     'targetL','targetR','blank'};
% trigNames = {'fastL-attL','fastL-attR','fastR-attL','fastR-attR','blank'};
trigNames = {'targetL','targetR'};

saveData = 0;
saveFigs = 0;

% load data header for plotting topologies
load data/data_hdr.mat

%% Get the data
trigMean = [];
trigData = [];
for iChSet = 1:numel(channelSets)
    allChannels = channelSets{iChSet};
    channels = setdiff(allChannels,badChannels);
    
    [trigM, triggers, Fs, trigD, trigEvents] =  rd_getData(filename, trigChan, channels, tstart, tstop);
    trigMean = cat(2,trigMean,trigM);
    trigData = cat(2,trigData,trigD);
end

nSamples = size(trigMean,1);
nChannels = size(trigMean,2);
nTrigs = size(trigMean,3);

%% Save the data
if saveData
    save(savename);
end

%% Exclude trials manually rejected with ft
if excludeTrialsFt
    % load trials_rejected variable from ft manual rejection
    load([dataDir '/mat/trials_rejected_erf.mat'])
    
    includedTrials = logical(ones(size(trigData,3),1));
    includedTrials(trials_rejected) = 0;
    
    trigMean = [];
    for iTrig = 1:nTrigs
        trigger = triggers(iTrig);
        w = trigEvents(:,2)==trigger & includedTrials;
        trigMean(:,:,iTrig) = mean(trigData(:,:,w),3);
    end
    trigData = trigData(:,:,includedTrials);
    trigEvents = trigEvents(includedTrials,:);
    
    % update figDir
    figDir = [figDir '_ft'];
end

%% ANALYZE THE MEAN TIMESERIES FOR EACH TRIGGER TYPE
%% Baseline
% baselinePeriod = -500:0;
baselinePeriod = t;
inBaseline = ismember(t,baselinePeriod);
baselineDC = squeeze(mean(mean(trigMean(inBaseline,:,:),1),3));
baselineTSeries = repmat(baselineDC,[size(trigMean,1),1,size(trigMean,3)]);

trigMean0 = trigMean;
trigMean = trigMean-baselineTSeries;

%% FFT on mean time series for each trigger type
% do the fft for each channel
nfft = 2^nextpow2(nSamples); % Next power of 2 from length of y
Y = fft(trigMean,nfft)/nSamples;
f = Fs/2*linspace(0,1,nfft/2+1);
amps = 2*abs(Y(1:nfft/2+1,:,:));

%% Plot trial average and single-sided amplitude spectrum
fH = [];
for iTrig = 1:nTrigs
    fH(iTrig) = figure;
    % time
    subplot(2,1,1)
    plot(repmat(t',1,nChannels), trigMean(:,:,iTrig))
    xlabel('time (ms)')
    ylabel('amplitude')
    title(trigNames{iTrig})
    % frequency
    subplot(2,1,2)
    plot(repmat(f',1,nChannels), amps(:,:,iTrig))
    xlim([1 200])
    ylim([0 20])
    xlabel('Frequency (Hz)')
    ylabel('|Y(f)|')
end

if saveFigs
    rd_saveAllFigs(fH, trigNames, 'plot_tsFFT', figDir);
end

%% Get some time series peaks
times = [228 475];
timeWindow = 5; % +/- this window value
for iT = 1:numel(times)
    time = times(iT);
    inTimeRange = t<time+timeWindow & t>time-timeWindow;
    peakTMeans(iT,:,:) = mean(trigMean(inTimeRange,:,:),1);
end

%% Convert to 157 channels
for iT = 1:numel(times)
    timeToPlot = times(iT);
    %     timeToPlot = 228;
    timeIdx = find(times==timeToPlot);
    peakTM = squeeze(peakTMeans(timeIdx,:,:))';
    inds = setdiff(0:156,badChannels)+1;
    peakTMeans157 = to157chan(peakTM,inds,'zeros');
    
    %% Plot on mesh
    % all conditions separately
    fH = [];
    for iTrig = 1:nTrigs
        sensorData = peakTMeans157(iTrig,:);
        figure
        fH(iTrig) = ssm_plotOnMesh(sensorData, trigNames{iTrig}, [], data_hdr, '2d');
        set(gca,'CLim',[-150 150])
    end
    if saveFigs
        figPrefix = sprintf('map_amp%dms', timeToPlot);
        rd_saveAllFigs(fH,trigNames,figPrefix,figDir)
    end
end

%% Weighted ERF
w = load([dataDir '/mat/weights.mat']);

inds = setdiff(0:156,badChannels)+1;
if ~isempty(w.weightChannels)
    for iChannel = 1:numel(inds)
        inWeightSet(iChannel) = any(inds(iChannel)==w.weightChannels);
    end
else
    inWeightSet = true(ones(size(inds)));
end

wCondNames = {'fastL','fastR'};
for iFreq = 1:numel(w.ssvefFreqs)
    wFreq = w.ssvefFreqs(iFreq);
    
    fH = figure;
    for iCond = 1:numel(wCondNames)
        wCond = wCondNames{iCond};
        
        for iTrig = 1:nTrigs
            weights = squeeze(w.w(:,iFreq,iCond));
            weightedERF(:,iTrig) = trigMean(:,inWeightSet,iTrig)*weights; % trigMean or abs(trigMean)
        end
        
        subplot(2,1,iCond)
        plot(t, weightedERF)
        ylim([-40 40])
        xlabel('time (ms)')
        ylabel('amplitude') % amplitude or abs(amplitude)
        legend(trigNames)
        title(sprintf('%s, %d Hz weights', wCond, wFreq))
    end
    if saveFigs
        figName = sprintf('ERF_%dHzWeights', wFreq); % ERF or weightedAbsERF
        rd_saveAllFigs(fH, {figName}, 'plot', figDir);
    end
end

%% ANALYZE SINGLE TRIAL DATA
%% Baseline
% baselinePeriod = -500:0;
baselinePeriod = t;
inBaseline = ismember(t,baselinePeriod);
baselineDC = squeeze(mean(mean(trigData(inBaseline,:,:),1),3));
baselineTSeries = repmat(baselineDC,[size(trigData,1),1,size(trigData,3)]);

trigData0 = trigData;
trigData = trigData-baselineTSeries;

%% FFT on mean time series for each trigger type
% do the fft for each channel
nfft = 2^nextpow2(nSamples); % Next power of 2 from length of y
dataY = fft(trigData,nfft)/nSamples;
f = Fs/2*linspace(0,1,nfft/2+1);
dataAmps = 2*abs(dataY(1:nfft/2+1,:,:));

%% Plot the FFT for every trial from a sample channel
sampleChannel = 1;
figure
plot(f, squeeze(dataAmps(:,1,:)))
xlim([1 200])
ylim([0 40])
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')

%% Take the frequency average across trials for each channel, for each trigger type
for iTrig = 1:nTrigs
    trigger = triggers(iTrig);
    w = trigEvents(:,2)==trigger;
    dataAmpsMean(:,:,iTrig) = mean(dataAmps(:,:,w),3);
end

%% Plot single-sided amplitude spectrum for single trial data average
% figure
fH = [];
for iTrig = 1:nTrigs
    fH(iTrig) = figure;
    % frequency
    plot(repmat(f',1,nChannels), dataAmpsMean(:,:,iTrig))
    xlim([1 200])
    ylim([0 40])
    xlabel('Frequency (Hz)')
    ylabel('|Y(f)|')
    title(sprintf('single trial FFT average, %s',trigNames{iTrig}))
end
if saveFigs
    rd_saveAllFigs(fH,trigNames,'plot_singleTrialFFT',figDir)
end

%% Get the component peaks
ssvefFreqs = [15 20 30 40];
nonSSVEFFreqs = [25 25 45]; % random freqs in a similar range
freqs = [ssvefFreqs nonSSVEFFreqs];
freqWindow = 0.2; % +/- this window value
for iF = 1:numel(freqs)
    freq = freqs(iF);
    inFreqRange = f<freq+freqWindow & f>freq-freqWindow;
    peakFreqs{iF} = f(inFreqRange);
    peakVals{iF} = dataAmpsMean(inFreqRange,:,:);
    peakMeans(iF,:,:) = mean(peakVals{iF},1);
end

%% Find the channels with high SSVEF SNR
peakSignal = peakMeans(1:numel(ssvefFreqs),:,:); % freqs x channels x triggers
peakNoise = mean(peakMeans(numel(ssvefFreqs)+1:end,:,:),1); % 1 x channels x triggers
peakSNR = (peakSignal./repmat(peakNoise,numel(ssvefFreqs),1));
peakSNRAll = mean(mean(peakSNR,3),1);

snrThresh = prctile(peakSNRAll,90);
highSNRChannels0 = peakSNRAll>snrThresh;

% for each flicker frequency, across trigger types
figure
imagesc(mean(peakSNR,3)')
set(gca,'XTick',1:numel(ssvefFreqs))
set(gca,'XTickLabel',ssvefFreqs)
xlabel('frequency (Hz)')
ylabel('channel')
title('SNR (ssvef/non-ssvef amplitude)')
colorbar
if saveFigs
    rd_saveAllFigs(gcf,{'singleTrialSSVEFSNR'},'im',figDir)
end

% across flicker frequencies
figure
hist(peakSNRAll)
xlabel('SNR (ssvef/non-ssvef amplitude)')
ylabel('number of channels')
if saveFigs
    rd_saveAllFigs(gcf,{'singleTrialSSVEFSNR'},'hist',figDir)
end

%% Convert SNR to 157 channels and plot on mesh
peakSNRAll157 = to157chan(peakSNRAll,inds,'zeros');

figure
fH = ssm_plotOnMesh(peakSNRAll157, 'peak SNR all flickers', [], data_hdr, '2d');
set(gca,'CLim',[0 4])
if saveFigs
    rd_saveAllFigs(gcf,{'singleTrialSSVEFSNR'},'map',figDir)
end

%% Plot ERF of high SNR channels
for iChSet = 1:numel(pickedChannels)
    highSNRChannels = pickedChannels{iChSet};
    fH = [];
    for iTrig = 1:nTrigs
        fH(iTrig) = figure;
        % time
        subplot(2,1,1)
        plot(repmat(t',1,nnz(highSNRChannels)), trigMean(:,highSNRChannels,iTrig))
        xlabel('time (ms)')
        ylabel('amplitude')
        title(trigNames{iTrig})
        % frequency
        subplot(2,1,2)
        plot(repmat(f',1,nnz(highSNRChannels)), dataAmpsMean(:,highSNRChannels,iTrig))
        xlim([1 200])
        ylim([0 40])
        xlabel('Frequency (Hz)')
        ylabel('|Y(f)|')
    end
    if saveFigs
        %     rd_saveAllFigs(fH,trigNames,'plot_tsFFTHighSNRChannels',figDir)
        name = sprintf('plot_tsFFTPickedChannels%s', pickedChannelsLabels{iChSet});
        rd_saveAllFigs(fH,trigNames,name,figDir)
    end
    
    % time, mean across channels
    figure
    plot(t',squeeze(mean(trigMean(:,highSNRChannels,:),2)))
    xlabel('time (ms)')
    ylabel('amplitude')
    ylim([-200 200])
    legend(trigNames)
    if saveFigs
        %     rd_saveAllFigs(gcf,{'erfHighSNRChannels'},'plot',figDir)
        name = sprintf('tsFFTPickedChannels%s', pickedChannelsLabels{iChSet});
        rd_saveAllFigs(gcf,{name},'plot',figDir)
    end
end

%% FFT of trigMean for high SNR channels
channels = highSNRChannels;
targetERF = squeeze(mean(trigMean(:,channels,:),2));

targetY = fft(targetERF,nfft)/nSamples;
targetAmps = 2*abs(targetY(1:nfft/2+1,:));

taper          = 'hanning';
foi            = 1:100;
t_ftimwin      = 10 ./ foi;
toi            = tstart/1000:0.01:tstop/1000;
tfAmps = [];
for iTrig = 1:nTrigs
    data = trigMean(:,channels,iTrig)'; % channels by samples
    [spectrum,ntaper,freqoi,timeoi] = ft_specest_mtmconvol(data, t/1000, ...
        'timeoi', toi, 'freqoi', foi, 'timwin', t_ftimwin, ...
        'taper', taper, 'dimord', 'chan_time_freqtap');
    specAmp = squeeze(mean(abs(spectrum),1)); % mean across channels
    tfAmps(:,:,iTrig) = specAmp';
end

% this is not really the right baseline. should compare to blank or target
% absent. however, this is conservative.
tfAmpsLRDiff = tfAmps(:,:,1) - tfAmps(:,:,2);

ssvefFreq = 30;
noiseDist = 5;
freqIdx = find(foi==ssvefFreq);
signal = tfAmps(freqIdx,:,:);
noise = tfAmps(freqIdx + [-noiseDist noiseDist],:,:);

% figures
ytick = 10:10:numel(foi);
xtick = 51:50:numel(toi);
clims = [0 15];
diffClims = [-10 10];
eventTimes = 0;

fH = [];
fH(1) = figure;
plot(f, targetAmps)
xlim([0 150])
xlabel('frequency (Hz)')
ylabel('amplitude')
legend(trigNames)
title(['channel' sprintf(' %d', channels)]);

fH(2) = figure;
for iTrig = 1:nTrigs
    subplot(1,2,iTrig)
    imagesc(tfAmps(:,:,iTrig),clims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    if iTrig==nTrigs
        xlabel('time (ms)')
        ylabel('frequency (Hz)')
    end
    title(trigNames{iTrig})
end
rd_supertitle(['channel' sprintf(' %d', channels)]);
rd_raiseAxis(gca);

fH(3) = figure;
imagesc(tfAmpsLRDiff,diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (ms)')
ylabel([trigNames{1} ' - ' trigNames{2}])
title(['channel' sprintf(' %d', channels)]);
colorbar

fH(4) = figure;
for iTrig = 1:nTrigs
    subplot(1,2,iTrig)
    hold on
    plot(toi, signal(:,:,iTrig),'k', 'LineWidth', 2)
    plot(toi, noise(:,:,iTrig))
    xlabel('time (ms)')
    ylabel('mtmconvol amplitude')
    title(trigNames{iTrig})
end
legend(num2str(ssvefFreq), num2str(ssvefFreq-noiseDist), num2str(ssvefFreq+noiseDist))


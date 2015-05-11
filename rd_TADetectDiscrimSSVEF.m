% rd_SSVEFAnalysis.m

%% Setup
% exptDir = '/Local/Users/denison/Data/TAPilot/MEG';
exptDir = '/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG';
sessionDir = 'R0817_20150504';
fileBase = 'R0817_TADeDi_5.4.15';
analStr = 'ebi'; % '', 'eti', etc.
excludeTrialsFt = 1;

dataDir = sprintf('%s/%s', exptDir, sessionDir);

switch analStr
    case ''
        filename = sprintf('%s/%s.sqd', dataDir, fileBase);
        savename = sprintf('%s/mat/%s_ssvef_workspace.mat', dataDir, fileBase);
        figDir = sprintf('%s/figures/raw', dataDir);
    otherwise
        filename = sprintf('%s/%s_%s.sqd', dataDir, fileBase, analStr);
        savename = sprintf('%s/mat/%s_%s_ssvef_workspace.mat', dataDir, fileBase, analStr);
        figDir = sprintf('%s/figures/%s', dataDir, analStr);
end
if ~exist(figDir,'dir')
    mkdir(figDir)
end

behavDir = sprintf('%s/Behavior/%s/analysis', exptDir(1:end-4), sessionDir);
behavFile = dir(sprintf('%s/*.mat', behavDir));
behav = load(sprintf('%s/%s', behavDir, behavFile.name));

% trigChan = 160:167;
trigChan = [160:163 166]; % stim/blank blocks
megChannels = 0:156;
channelSets = {0:39,40:79,80:119,120:156};

switch sessionDir
    case 'R0890_20140806'
        % badChannels = [10 11 115]; % R0890, 48-->49, 150-->152
        badChannels = [10 11 115 49 152]; % R0890
    case 'R0817_20140820'
        badChannels = [115 152]; % R0817
        weightChannels = sort(unique([59 92 10 60 15 14 32 2 51 1 50 39 7 24 55 103 98 8]));
    otherwise
        badChannels = [];
        fprintf('sessionDir not found ... no bad channels selected\n\n')
end
% badChannels = [];

Fs = 1000;
tstart = -500; % ms % TAPilot: 1000
tstop = 4400; % ms % TAPilot: 6500
t = tstart:tstop;

% trigNames = {'fastL-attL','fastL-attR','fastR-attL','fastR-attR',...
%     'targetL','targetR','blank'};
% trigNames = {'fastL-attL','fastL-attR','fastR-attL','fastR-attR','blank'};
trigNames = {'attT1-T1p-T2p','attT2-T1p-T2p','attT1-T1a-T2p','attT2-T1a-T2p',...
    'attT1-T1p-T2a','attT2-T1p-T2a','attT1-T1a-T2a','attT2-T1a-T2a','blank'};

saveData = 0;
saveFigs = 0;

% load data header for plotting topologies
load data/data_hdr.mat

%% Get the data
trigData = [];
% trigMean = [];
% triggers = [];
for iChSet = 1:numel(channelSets)
    allChannels = channelSets{iChSet};
    channels = setdiff(allChannels,badChannels);
    
    [trigM, triggers, Fs, trigD, trigEvents] =  rd_getData(filename, trigChan, channels, tstart, tstop);
%     trigMean = cat(2,trigMean,trigM);
    trigData = cat(2,trigData,trigD);
end
nSamples = size(trigMean,1);
nChannels = size(trigMean,2);

%% Save the data
if saveData
    save(savename);
end

%% Baseline
% baselinePeriod = -500:0;
baselinePeriod = t;
inBaseline = ismember(t,baselinePeriod);
baselineDC = mean(trigData(inBaseline,:,:),1);
baselineTSeries = repmat(baselineDC,[size(trigData,1),1,1]);

trigData0 = trigData;
trigData = trigData-baselineTSeries;

%% Organize trials into conditions
cueCondIdx = strcmp(behav.responseData_labels, 'cue condition');
t1CondIdx = strcmp(behav.responseData_labels, 'target type T1');
t2CondIdx = strcmp(behav.responseData_labels, 'target type T2');

blankCond = 1;
cueConds = {[2 3], [4 5]}; % cue T1, cue T2
t1Conds = {[1 2], 0}; % present, absent
t2Conds = {[1 2], 0}; % present, absent

condData = [];
for iCue = 1:numel(cueConds)
    vals = cueConds{iCue};
    wCue = [];
    for iEl = 1:numel(vals)
        wCue(:,iEl) = behav.responseData_all(:,cueCondIdx) == vals(iEl);
    end
    for iT1 = 1:numel(t1Conds)
        vals = t1Conds{iT1};
        wT1 = [];
        for iEl = 1:numel(vals)
            wT1(:,iEl) = behav.responseData_all(:,t1CondIdx) == vals(iEl);
        end
        for iT2 = 1:numel(t2Conds)
            vals = t2Conds{iT2};
            wT2 = [];
            for iEl = 1:numel(vals)
                wT2(:,iEl) = behav.responseData_all(:,t2CondIdx) == vals(iEl);
            end
            
            w = sum([wCue wT1 wT2],2)==3;
            condData(:,:,:,iCue,iT1,iT2) = trigData(:,:,w);
        end
    end
end

wBlank = behav.responseData_all(:,cueCondIdx) == blankCond;
blankData = trigData(:,:,wBlank);

% mean across trials
condDataMean = squeeze(mean(condData,3));
blankDataMean = squeeze(mean(blankData,3));

% let trigMean have the conditions 1-9 in the third dimension
trigMean = condDataMean(:,:,:);
trigMean(:,:,end+1) = blankDataMean;
nTrigs = size(trigMean,3);

%% Exclude trials manually rejected with ft
% if excludeTrialsFt
%     % load trials_rejected variable from ft manual rejection
%     load([dataDir '/mat/trials_rejected_ssvef.mat'])
%     
%     includedTrials = logical(ones(size(trigData,3),1));
%     includedTrials(trials_rejected) = 0;
%     
%     trigMean = [];
%     for iTrig = 1:nTrigs
%         trigger = triggers(iTrig);
%         w = trigEvents(:,2)==trigger & includedTrials;
%         trigMean(:,:,iTrig) = mean(trigData(:,:,w),3);
%     end
%     trigData = trigData(:,:,includedTrials);
%     trigEvents = trigEvents(includedTrials,:);
%     
%     % update figDir
%     figDir = [figDir '_ft'];
% end

%% Find noisy channels
varCutoff = 100;
trigVar = std(trigMean(:,:,1));
figure
hist(trigVar)
noisyChannels = trigVar>varCutoff;

%% FFT on mean time series for each trigger type
% do the fft for each channel
nfft = 2^nextpow2(nSamples); % Next power of 2 from length of y
Y = fft(trigMean,nfft)/nSamples; % Scale by number of samples
f = Fs/2*linspace(0,1,nfft/2+1); % Fs/2 is the maximum frequency that can be measured
amps = 2*abs(Y(1:nfft/2+1,:,:)); % Multiply by 2 since only half the energy is in the positive half of the spectrum?

%% Plot trial average and single-sided amplitude spectrum
% figure
fH = [];
for iTrig = 1:nTrigs
    fH(iTrig) = figure;
    % time
    subplot(2,1,1)
    %     hold on
    plot(repmat(t',1,nChannels), trigMean(:,:,iTrig))
    xlabel('time (ms)')
    ylabel('amplitude')
    title(trigNames{iTrig})
    % frequency
    subplot(2,1,2)
    %     hold on
    plot(repmat(f',1,nChannels), amps(:,:,iTrig))
    xlim([1 200])
    ylim([0 20])
    xlabel('Frequency (Hz)')
    ylabel('|Y(f)|')
end

if saveFigs
    rd_saveAllFigs(fH, trigNames, 'plot_tsFFT', figDir);
end

%% Get the component peaks
ssvefFreqs = [15 20 30 40];
freqWindow = 0.2; % +/- this window value
for iF = 1:numel(ssvefFreqs)
    freq = ssvefFreqs(iF);
    inFreqRange = f<freq+freqWindow & f>freq-freqWindow;
    peakFreqs{iF} = f(inFreqRange);
    peakVals{iF} = amps(inFreqRange,:,:);
    peakMeans(iF,:,:) = mean(peakVals{iF},1);
end

%% Plot peak freq image
fH = [];
for iF = 1:numel(ssvefFreqs)
    freq = ssvefFreqs(iF);
    fH(iF) = figure;
    imagesc(squeeze(peakMeans(iF,:,:)))
    title(sprintf('frequency = %d',freq))
    freqNames{iF} = sprintf('peakAmp%dHz', freq);
end
if saveFigs
    rd_saveAllFigs(fH, freqNames, 'im', figDir)
end

%% Convert to 157 channels
for iF = 1:numel(ssvefFreqs)
    freqToPlot = ssvefFreqs(iF);
%     freqToPlot = 40;
    freqIdx = find(ssvefFreqs==freqToPlot);
    peakM = squeeze(peakMeans(freqIdx,:,:))';
    inds = setdiff(0:156,badChannels)+1;
    peakMeans157 = to157chan(peakM,inds,'zeros');
    
    stimAve = [ones(1, nTrigs-1) 0]./(nTrigs-1);
    peakStimAve157 = (peakMeans157'*stimAve')';

    %% Plot on mesh
    % all conditions separately
    fH = [];
    for iTrig = 1:nTrigs
        sensorData = peakMeans157(iTrig,:);
        figure
        fH(iTrig) = ssm_plotOnMesh(sensorData, trigNames{iTrig}, [], data_hdr, '2d');
        set(gca,'CLim',[0 15])
    end
    
    % stim average
    figure
    fH(end+1) = ssm_plotOnMesh(peakStimAve157, 'stim average', [], data_hdr, '2d');
    set(gca,'CLim',[0 15])

    % save figs
    if saveFigs
        figNames = [trigNames {'LRDiff','AttInOutDiff','LStimAttInOutDiff','RStimAttInOutDiff','AttEffectLRDiff'}];
        figPrefix = sprintf('map_ssvef%dHz', freqToPlot);
        rd_saveAllFigs(fH,figNames,figPrefix,figDir)
    end
end

%% Find the channels with high SSVEF SNR
peakSignal = mean(peakMeans(:,:,1:4),3); % freqs x channels
peakNoise = mean(peakMeans(:,:,5),3); % freqs x channels
peakSNR = (peakSignal./peakNoise)';
peakSNRAllFlickers = mean(peakSNR,2);

figure
imagesc(peakSNR)

figure
hist(peakSNRAllFlickers)
xlabel('peak SNR (mean across flicker frequencies)')
ylabel('number of channels')

if saveFigs
    rd_saveAllFigs(gcf, {'peakSNRAllFlickers'}, 'hist', figDir)
end

%% Plot selected channels for all ssvef freqs, all conditions
% channelsToPlot = find(peakSNRAllFlickers>10);
channelsToPlot = [25 13];
for iF = 1:numel(ssvefFreqs)
    freq = ssvefFreqs(iF);
    figure
    bar(squeeze(peakMeans(iF,channelsToPlot,:)))
    set(gca,'XTickLabel',channelsToPlot)
    legend(trigNames)
    title(sprintf('%d Hz peak', freq))
end

%% Plotting setup
plotOrder = [1 5 3 7 2 6 4 8 9];
extendedMap = flipud(lbmap(nTrigs-1+4,'RedBlue'));
selectedMap = extendedMap([1:(nTrigs-1)/2 (end-(nTrigs-1)/2)+1:end],:);
trigColors = [selectedMap; 0 0 0];
trigBlue = mean(selectedMap(1:(nTrigs-1)/2,:));
trigRed = mean(selectedMap((end-(nTrigs-1)/2)+1:end,:));

eventTimes = [0 1000 1600 2600];

%% Wavelet on average across trials
channels = 25; % [13 25]
ssvefFreq = 30;
for iTrig = 1:nTrigs
    data = trigMean(:,channels,iTrig)'; % channels by samples
    [spectrum,freqoi,timeoi] = ft_specest_wavelet(data, t/1000);
    specAmp = abs(squeeze(spectrum));
    
    freqIdx = find(abs(freqoi-ssvefFreq) == min((abs(freqoi-ssvefFreq))));
    
    wAmp = squeeze(specAmp(:,freqIdx,:));
    wAmpNorm = wAmp./nanmean(nanmean(wAmp(:,1:500)));
    wAmps0(:,:,iTrig) = wAmpNorm';
end
wAmps = squeeze(mean(wAmps0,2)); % mean across channels

figure
set(gca,'ColorOrder',trigColors)
hold all
plot(t, wAmps(:,plotOrder))
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
plot(t, mean(wAmps(:,plotOrder(1:(nTrigs-1)/2)),2),'color',trigBlue,'LineWidth',4)
plot(t, mean(wAmps(:,plotOrder(end-(nTrigs-1)/2):end-1),2),'color',trigRed,'LineWidth',4)
legend(trigNames(plotOrder))
xlabel('time (ms)')
ylabel('wavelet amp')
title(['channels' sprintf(' %d', channels)])

%% Hilbert on average across trials
channel = 25;
ssvefFreq = 30;
Fbp = ssvefFreq + [-1.6 1.6];
for iTrig = 1:nTrigs
    data = trigMean(:,channel,iTrig)'; % channels by samples
    dataF = ft_preproc_bandpassfilter(data,Fs,Fbp);
    dataFH = abs(hilbert(dataF));
    hAmps(:,iTrig) = dataFH;
end

figure
set(gca,'ColorOrder',trigColors)
hold all
plot(t, hAmps(:,plotOrder))
legend(trigNames(plotOrder))
title(['channel' sprintf(' %d', channel)])
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end

%% Filter single trial data
ssvefFreq = 30;
Fbp = ssvefFreq + [-1.6 1.6];
fprintf('\nFiltering single trial data ...\n')
for iTrial = 1:size(condData,3)
    fprintf('trial %d\n', iTrial)
    for iCue = 1:numel(cueConds)
        for iT1 = 1:numel(t1Conds)
            for iT2 = 1:numel(t2Conds)
                data = condData(:,:,iTrial,iCue,iT1,iT2)'; % channels by samples
                dataF = ft_preproc_bandpassfilter(data,Fs,Fbp);
                condDataF(:,:,iTrial,iCue,iT1,iT2) = dataF';
            end
        end
    end
end

for iTrial = 1:size(blankData,3)
    data = blankData(:,:,iTrial)';
    dataF = ft_preproc_bandpassfilter(data,Fs,Fbp);
    blankDataF(:,:,iTrial) = dataF';
end

% mean across trials
condDataFMean = squeeze(mean(condDataF,3));
blankDataFMean = squeeze(mean(blankDataF,3));

% let trigFMean have the conditions 1-9 in the third dimension
trigFMean = condDataFMean(:,:,:);
trigFMean(:,:,end+1) = blankDataFMean;

%% Hilbert on filtered single trial means
for iChan = 1:nChannels
    for iTrig = 1:nTrigs
        dataF = trigFMean(:,iChan,iTrig);
        dataFH = abs(hilbert(dataF));
        hTrigFMean(:,iChan,iTrig) = dataFH;
    end
end

channel = 13;
figure
set(gca,'ColorOrder',trigColors)
hold all
plot(t, squeeze(hTrigFMean(:,channel,plotOrder)))
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
plot(t, mean(squeeze(hTrigFMean(:,channel,plotOrder(1:(nTrigs-1)/2))),2),'color',trigBlue,'LineWidth',4)
plot(t, mean(squeeze(hTrigFMean(:,channel,plotOrder(end-(nTrigs-1)/2):end-1)),2),'color',trigRed,'LineWidth',4)
legend(trigNames(plotOrder))
xlabel('time (ms)')
ylabel('hilbert amp')


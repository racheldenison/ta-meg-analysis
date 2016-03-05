function rd_TADetectDiscrimSSVEF1(sessionDir)
% rd_TADetectDiscrimSSVEF1.m

%% Setup
exptDir = '/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG';
% sessionDir = 'R0817_20150504';
fileBase = sessionDirToFileBase(sessionDir);
analStr = 'ebi'; % '', 'eti', 'ebi', etc.
excludeTrialsFt = 1;
excludeSaturatedEpochs = 0;

dataDir = sprintf('%s/%s', exptDir, sessionDir);
matDir = sprintf('%s/mat', dataDir);

switch analStr
    case ''
        filename = sprintf('%s/%s.sqd', dataDir, fileBase);
        savename = sprintf('%s/%s_ssvef_workspace.mat', matDir, fileBase);
        figDir = sprintf('%s/figures/raw', dataDir);
    otherwise
        filename = sprintf('%s/%s_%s.sqd', dataDir, fileBase, analStr);
        savename = sprintf('%s/%s_%s_ssvef_workspace.mat', matDir, fileBase, analStr);
        figDir = sprintf('%s/figures/%s', dataDir, analStr);
end

behavDir = sprintf('%s/Behavior/%s/analysis', exptDir(1:end-4), sessionDir);
behavFile = dir(sprintf('%s/*.mat', behavDir));
behav = load(sprintf('%s/%s', behavDir, behavFile.name));

% trigChan = 160:167;
trigChan = [160:163 166]; % stim/blank blocks
megChannels = 0:156;
channelSets = {0:39,40:79,80:119,120:156};
% for checking triggers:
% megChannels = 160:167;
% channelSets = {megChannels};

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

Fs = 1000;
tstart = -500; % ms 
tstop = 3600; % ms 
t = tstart:tstop;

eventTimes = [0 500 1500 2100 3100];

trigNames = {'attT1-T1p-T2p','attT2-T1p-T2p','attT1-T1a-T2p','attT2-T1a-T2p',...
    'attT1-T1p-T2a','attT2-T1p-T2a','attT1-T1a-T2a','attT2-T1a-T2a','blank'};
% for checking triggers:
% tn = {'1-1','1-2','2-1','2-2','abs','pres','blank','cue'};

saveData = 0;
saveFigs = 0;

% load data header for plotting topologies
load data/data_hdr.mat

%% Get the data
if ~exist(savename, 'file')
    trigData = [];
    for iChSet = 1:numel(channelSets)
        allChannels = channelSets{iChSet};
        channels = setdiff(allChannels,badChannels);
        
        [trigM, triggers, Fs, trigD, trigEvents] =  rd_getData(filename, trigChan, channels, tstart, tstop);
        trigData = cat(2,trigData,trigD);
    end
    nSamples = size(trigData,1);
    nChannels = size(trigData,2);
else
    load(savename)
end

%% Save the data
if saveData
    if ~exist(matDir,'dir')
        mkdir(matDir)
        
        prepFile = sprintf('%s/prep/trials_rejected.mat',dataDir);
        if exist(prepFile,'file')
            movefile(prepFile,matDir)
        end 
    end
    save(savename, '-v7.3');
end

%% Find saturated channels and trials in raw data
if strcmp(analStr, '')
    saturatedChannelEpochs = rd_findSaturatedChannelEpochs(trigData);
    if size(saturatedChannelEpochs, 2)~=41*14
        fprintf('\nMake sure we are taking the right trials!\n')
        saturatedChannelEpochs = saturatedChannelEpochs(:,42:end);
        fprintf('\nNew size: [%d %d]\n\n', size(saturatedChannelEpochs))
    end
    save(sprintf('%s/saturated_channel_epochs.mat', matDir), 'saturatedChannelEpochs');
    if saveFigs
        rd_saveAllFigs(gcf, {'saturatedChannelEpochs'}, 'im', figDir);
    end
end

%% Baseline
baselinePeriod = t;
inBaseline = ismember(t,baselinePeriod);
baselineDC = mean(trigData(inBaseline,:,:),1);
baselineTSeries = repmat(baselineDC,[size(trigData,1),1,1]);

trigData = trigData-baselineTSeries;

%% Excluded saturated channel epochs
if excludeSaturatedEpochs
    load([matDir '/saturated_channel_epochs.mat'])
    trigData(:,saturatedChannelEpochs) = NaN;
end

%% Exclude trials manually rejected with ft
if excludeTrialsFt
    % load trials_rejected variable from ft manual rejection
    load([matDir '/trials_rejected.mat'])

    trigData(:,:,trials_rejected) = NaN;
    
    % update figDir
    figDir = [figDir '_ft'];
end

%% Fig dir
if ~exist(figDir,'dir')
    mkdir(figDir)
end

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
condDataMean = squeeze(nanmean(condData,3));
blankDataMean = squeeze(nanmean(blankData,3));

% let trigMean have the conditions 1-9 in the third dimension
trigMean = condDataMean(:,:,:);
trigMean(:,:,end+1) = blankDataMean;
nTrigs = size(trigMean,3);

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
    hold on
    plot(repmat(t',1,nChannels), trigMean(:,:,iTrig))
    for iEv = 1:numel(eventTimes)
        vline(eventTimes(iEv),'k');
    end
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

% Average across conditions
figure
% time
subplot(2,1,1)
hold on
plot(repmat(t',1,nChannels), mean(trigMean,3))
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
xlabel('time (ms)')
ylabel('amplitude')
title('stim average')
% frequency
subplot(2,1,2)
%     hold on
plot(repmat(f',1,nChannels), mean(amps,3))
xlim([1 200])
ylim([0 20])
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')

if saveFigs
    rd_saveAllFigs(gcf, {'stimAve'}, 'plot_tsFFT', figDir)
end

%% Get the component peaks
ssvefFreqs = [15 20 30 40];
% ssvefFreqs = [10 50 57.62 60];
% ssvefFreqs = [1 2.68 4 5 6.7 8.4 13.43];
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

%% Plot peaks for stim ave, save channelsRanked
for ssvefFreq = [30 40]
    % ssvefFreq = 40;
    peakMeansStimAve = squeeze(mean(peakMeans(ssvefFreqs==ssvefFreq,:,1:end-1),3));
    peakMeansBlank = squeeze(mean(peakMeans(ssvefFreqs==ssvefFreq,:,end),3));
    [channelsRankedAmps, channelsRanked] = sort(peakMeansStimAve,2,'descend');
    channelsRanked(isnan(channelsRankedAmps)) = [];
    channelsRankedAmps(isnan(channelsRankedAmps)) = [];

    figure
    bar(peakMeansStimAve)
    text(120, channelsRankedAmps(1)-1, sprintf('top 5 channels:\n%s', num2str(channelsRanked(1:5))))
    xlabel('channel')
    ylabel('SSVEF amplitude')
    title(sprintf('%d Hz', ssvefFreq))
    
    if saveData
        if strcmp(analStr,'')
            channelsFileName = sprintf('%s/channels_%dHz.mat', matDir, ssvefFreq);
        else
            channelsFileName = sprintf('%s/channels_%dHz_%s.mat', matDir, ssvefFreq, analStr);
        end
        save(channelsFileName,'channelsRanked','channelsRankedAmps','peakMeansStimAve','peakMeansBlank')
    end
    if saveFigs
        figPrefix = 'bar';
        rd_saveAllFigs(gcf, {sprintf('channelStimAveAmp_%dHz', ssvefFreq)}, figPrefix, figDir)
    end
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
        set(gca,'CLim',[0 20])
    end
    
    % stim average
    figure
    fH(end+1) = ssm_plotOnMesh(peakStimAve157, 'stim average', [], data_hdr, '2d');
    set(gca,'CLim',[0 20])

    % save figs
    if saveFigs
        figNames = [trigNames {'StimAve'}];
        figPrefix = sprintf('map_ssvef%dHz', freqToPlot);
        rd_saveAllFigs(fH,figNames,figPrefix,figDir)
    end
end


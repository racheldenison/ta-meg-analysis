% rd_TADetectDiscrimSSVEF.m

%% Setup
% exptDir = '/Local/Users/denison/Data/TAPilot/MEG';
exptDir = '/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG';
sessionDir = 'R0988_20150904';
fileBase = 'R0988_TADeDi_r1-8_9.4.15';
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
% badChannels = [];

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
    % trigMean = [];
    % triggers = [];
    for iChSet = 1:numel(channelSets)
        allChannels = channelSets{iChSet};
        channels = setdiff(allChannels,badChannels);
        
        [trigM, triggers, Fs, trigD, trigEvents] =  rd_getData(filename, trigChan, channels, tstart, tstop);
        %     trigMean = cat(2,trigMean,trigM);
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
% baselinePeriod = -500:0;
baselinePeriod = t;
inBaseline = ismember(t,baselinePeriod);
baselineDC = mean(trigData(inBaseline,:,:),1);
baselineTSeries = repmat(baselineDC,[size(trigData,1),1,1]);

% trigData0 = trigData;
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
    
%     includedTrials = true(size(trigData,3),1);
%     includedTrials(trials_rejected) = 0;
    
%     trigMean = [];
%     for iTrig = 1:nTrigs
%         trigger = triggers(iTrig);
%         w = trigEvents(:,2)==trigger & includedTrials;
%         trigMean(:,:,iTrig) = mean(trigData(:,:,w),3);
%     end
%     trigData = trigData(:,:,includedTrials);
%     trigEvents = trigEvents(includedTrials,:);
    trigData(:,:,trials_rejected) = NaN;
    
    % update figDir
    figDir = [figDir '_ft'];
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
ssvefFreq = 40;
peakMeansStimAve = squeeze(mean(peakMeans(ssvefFreqs==ssvefFreq,:,1:end-1),3));
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
    save(channelsFileName,'channelsRanked','channelsRankedAmps')
end
if saveFigs
    figPrefix = 'bar';
    rd_saveAllFigs(gcf, {sprintf('channelStimAveAmp_%dHz', ssvefFreq)}, figPrefix, figDir)
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

%% Plotting setup
plotOrder = [1 5 3 7 2 6 4 8 9];
extendedMap = flipud(lbmap(nTrigs-1+4,'RedBlue'));
selectedMap = extendedMap([1:(nTrigs-1)/2 (end-(nTrigs-1)/2)+1:end],:);
trigColors = [selectedMap; 0 0 0];
trigBlue = mean(selectedMap(1:(nTrigs-1)/2,:));
trigRed = mean(selectedMap((end-(nTrigs-1)/2)+1:end,:));

% trigColorsPA4 = [107 76 154; 62 150 81; 57 106 177; 218 124 48]./255;
trigColorsPA4 = [.52 .37 .75; .31 .74 .40; .27 .51 .84; 1.0 .57 .22];

tsFigPos = [0 500 1250 375];
ts2FigPos = [0 500 1100 600];
ts3FigPos = [0 500 1100 900];
condFigPos = [250 300 750 650];
tf9FigPos = [0 250 1280 580];
tf3FigPos = [200 475 1000 275];

set(0,'defaultLineLineWidth',1)

%% Time series and FFT for single channel
channel = 60;
figure
set(gcf,'Position',ts2FigPos)

% time
subplot(3,1,1)
set(gca,'ColorOrder',trigColors)
hold all
plot(t, squeeze(trigMean(:,channel,plotOrder)))
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
% plot(t, squeeze(mean(trigMean(:,channel,plotOrder(1:(nTrigs-1)/2)),3)),'color',trigBlue,'LineWidth',4)
% plot(t, squeeze(mean(trigMean(:,channel,plotOrder(end-(nTrigs-1)/2):end-1),3)),'color',trigRed,'LineWidth',4)
xlim([t(1) t(end)])
xlabel('time (ms)')
ylabel('amplitude')
title(sprintf('channel %d', channel))

% frequency
subplot(3,1,2)
set(gca,'ColorOrder',trigColors)
hold all
plot(repmat(f',1,nTrigs), squeeze(amps(:,channel,plotOrder)))
xlim([1 200])
ylim([0 20])
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
legend(trigNames(plotOrder))

subplot(3,1,3)
set(gca,'ColorOrder',[.66 .5 .78; trigColors(end,:)])
hold all
plot(f, squeeze(mean(amps(:,channel,1:end-1),3)))
plot(f, squeeze(amps(:,channel,end)))
xlim([1 200])
ylim([0 20])
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
legend('stim average','blank')

if saveFigs
    figPrefix = sprintf('plot_ch%d', channel);
    rd_saveAllFigs(gcf, {'tsFFT'}, figPrefix, figDir)
end

%% Trial average for target present vs. absent, for a single channel
pp = mean(trigMean(:,:,1:2),3);
pa = mean(trigMean(:,:,5:6),3);
ap = mean(trigMean(:,:,3:4),3);
aa = mean(trigMean(:,:,7:8),3);

channel = 60;
targetWindow = [-100 500];
t1Window = t>=eventTimes(3) + targetWindow(1) & t<=eventTimes(3) + targetWindow(2);
t2Window = t>=eventTimes(4) + targetWindow(1) & t<=eventTimes(4) + targetWindow(2);
targetPA(1,:,1) = pp(t1Window, channel);
targetPA(2,:,1) = pp(t2Window, channel);
targetPA(3,:,1) = pa(t1Window, channel);
targetPA(4,:,1) = ap(t2Window, channel);
targetPA(1,:,2) = aa(t1Window, channel);
targetPA(2,:,2) = aa(t2Window, channel);
targetPA(3,:,2) = ap(t1Window, channel);
targetPA(4,:,2) = pa(t2Window, channel);

targetPADiff = targetPA(:,:,1)-targetPA(:,:,2);

targetNfft = 2^nextpow2(diff(targetWindow)+1);
targetY = fft(mean(targetPADiff),targetNfft)/(diff(targetWindow)+1);
targetF = Fs/2*linspace(0,1,targetNfft/2+1);
targetAmps = 2*abs(targetY(1:targetNfft/2+1));

names = {'target present','target absent'};
fH = [];
fH(1) = figure;
set(gcf,'Position',ts2FigPos)
for iPA = 1:2
    subplot(2,1,iPA)
    plot(targetWindow(1):targetWindow(2), targetPA(:,:,iPA))
    vline(0,'k');
    if iPA==1
        legend('xx1','xx2','xo1','ox2')
    end
    xlabel('time (ms)')
    ylabel('amplitude')
    title(names{iPA})
end
rd_supertitle(sprintf('channel %d', channel))
rd_raiseAxis(gca);

fH(2) = figure;
set(gcf,'Position',ts3FigPos)
subplot(3,1,1)
plot(targetWindow(1):targetWindow(2), targetPADiff)
xlabel('time (ms)')
ylabel('\Delta amplitude')
title('target present - absent')
subplot(3,1,2)
plot(targetWindow(1):targetWindow(2), mean(targetPADiff,1), 'k');
xlabel('time (ms)')
ylabel('\Delta amplitude')
subplot(3,1,3)
plot(targetF, targetAmps)
xlim([0 150])
xlabel('frequency (Hz)')
ylabel('\Delta amplitude')
rd_supertitle(sprintf('channel %d', channel))
rd_raiseAxis(gca);

if saveFigs
    figPrefix = sprintf('plot_ch%d', channel);
    rd_saveAllFigs(fH, {'targetPATrialAve','targetPATrialAveDiff'}, figPrefix, figDir)
end

%% Wavelet on average across trials
channels = 60; % [14 23 25 26 60], [13 14 23 25 43], [7 8 13 20 36]
ssvefFreq = 30;
width = 12; % 12 for 30 Hz, 16 for 40 Hz gives 127 ms duration, 5 Hz bandwidth
wBaselineWindow = [-500 0]; % [-300 -200];
wBaselineWindowIdx = find(t==wBaselineWindow(1)):find(t==wBaselineWindow(2));

% all frequencies
% wAmps0 = [];
% for iTrig = 1:nTrigs
%     data = trigMean(:,channels,iTrig)'; % channels by samples
%     [spectrum,freqoi,timeoi] = ft_specest_wavelet(data, t/1000, 'width', width);
%     specAmp = abs(squeeze(spectrum));
%     
%     freqIdx = find(abs(freqoi-ssvefFreq) == min((abs(freqoi-ssvefFreq))));
%     
%     if numel(size(specAmp))==3 % if three-dimensional
%         wAmp = squeeze(specAmp(:,freqIdx,:));
%     else
%         wAmp = squeeze(specAmp(freqIdx,:));
%     end
%     wAmpNorm = wAmp./nanmean(nanmean(wAmp(:,wBaselineWindowIdx)))-1;
%     wAmps0(:,:,iTrig) = wAmpNorm';
% end
% wAmps = squeeze(mean(wAmps0,2)); % mean across channels

% only frequency of interest
wAmps0 = [];
foi = ssvefFreq;
for iTrig = 1:nTrigs
    data = trigMean(:,channels,iTrig)'; % channels by samples
    [spectrum,freqoi,timeoi] = ft_specest_wavelet(data, t/1000, 'freqoi', foi, 'width', width);
    specAmp = abs(squeeze(spectrum));

    if all(size(specAmp)>1) % if two-dimensional
        wAmp = specAmp;
    else
        wAmp = specAmp';
    end
    wAmpNorm = wAmp./nanmean(nanmean(wAmp(:,wBaselineWindowIdx)))-1;
    wAmps0(:,:,iTrig) = wAmpNorm';
end
wAmps = squeeze(mean(wAmps0,2)); % mean across channels

fH = [];
fH(1) = figure;
set(gcf,'Position',tsFigPos)
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
title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels)])

% condition subplots
fH(2) = figure;
set(gcf,'Position',condFigPos)
for iTrig = 1:(nTrigs-1)/2
    subplot((nTrigs-1)/2,1,iTrig)
    set(gca,'ColorOrder',[trigBlue; trigRed])
    hold all
    plot(t, wAmps(:,iTrig*2-1:iTrig*2))
    legend(trigNames{iTrig*2-1:iTrig*2})
    ylim([-1 2.5])
    for iEv = 1:numel(eventTimes)
        vline(eventTimes(iEv),'k');
    end
    if iTrig==1
        title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels)])
    end
end
xlabel('time (ms)')
ylabel('wavelet amp')

% present/absent
fH(3) = figure;
set(gcf,'Position',tsFigPos)
hold on
for iTrig = 1:(nTrigs-1)/2 
    p1 = plot(t, mean(wAmps(:,iTrig*2-1:iTrig*2),2));
    set(p1, 'Color', trigColorsPA4(iTrig,:), 'LineWidth', 1.5)
end
ylim([-1 2.5])
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
legend('T1p-T2p','T1a-T2p','T1p-T2a','T1a-T2a')
xlabel('time (ms)')
ylabel('wavelet amp')
title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels)])

if saveFigs
    figPrefix = ['plot_ch' sprintf('%d_', channels) sprintf('%dHz', ssvefFreq)];
    rd_saveAllFigs(fH, {'waveletTrialAve','waveletTrialAveByCond','waveletTrialAvePA'}, figPrefix, figDir)
end

%% Wavelet on single trials
channel = 60;
wAmpsCond0 = [];

% all frequencies
% for iCue = 1:numel(cueConds)
%     for iT1 = 1:numel(t1Conds)
%         for iT2 = 1:numel(t2Conds)
%             data = squeeze(condData(:,channel,:,iCue,iT1,iT2))'; % trials by samples
%             [spectrum,freqoi,timeoi] = ft_specest_wavelet(data, t/1000);
%             specAmp = abs(squeeze(spectrum));
%             
%             freqIdx = find(abs(freqoi-ssvefFreq) == min((abs(freqoi-ssvefFreq))));
%             
%             wAmp = squeeze(specAmp(:,freqIdx,:));
%             wAmpNorm = wAmp./nanmean(nanmean(wAmp(:,wBaselineWindowIdx)))-1;
%             wAmpsCond0(:,:,iCue,iT1,iT2) = wAmpNorm';
%         end
%     end
% end
% 
% % blank
% data = squeeze(blankData(:,channel,:))'; % trials by samples
% [spectrum,freqoi,timeoi] = ft_specest_wavelet(data, t/1000);
% specAmp = abs(squeeze(spectrum));
% freqIdx = find(abs(freqoi-ssvefFreq) == min((abs(freqoi-ssvefFreq))));
% wAmp = squeeze(specAmp(:,freqIdx,:));
% wAmpNorm = wAmp./nanmean(nanmean(wAmp(:,wBaselineWindowIdx)))-1;
% wAmpsBlank0 = wAmpNorm';

% only frequency of interest
foi = ssvefFreq;
for iCue = 1:numel(cueConds)
    for iT1 = 1:numel(t1Conds)
        for iT2 = 1:numel(t2Conds)
            data = squeeze(condData(:,channel,:,iCue,iT1,iT2))'; % trials by samples
            [spectrum,freqoi,timeoi] = ft_specest_wavelet(data, t/1000, 'freqoi', foi, 'width', width);
            specAmp = abs(squeeze(spectrum));

            wAmp = specAmp;
            wAmpNorm = wAmp./nanmean(nanmean(wAmp(:,wBaselineWindowIdx)))-1;
            wAmpsCond0(:,:,iCue,iT1,iT2) = wAmpNorm';
        end
    end
end

% blank
data = squeeze(blankData(:,channel,:))'; % trials by samples
[spectrum,freqoi,timeoi] = ft_specest_wavelet(data, t/1000, 'freqoi', foi, 'width', width);
specAmp = abs(squeeze(spectrum));
wAmp = specAmp;
wAmpNorm = wAmp./nanmean(nanmean(wAmp(:,wBaselineWindowIdx)))-1;
wAmpsBlank0 = wAmpNorm';

% mean across trials
wAmpsCond0Mean = squeeze(nanmean(wAmpsCond0,2));
wAmpsBlank0Mean = squeeze(nanmean(wAmpsBlank0,2));

% let trigFMean have the conditions 1-9 in the third dimension
wAmpsTrigMean = wAmpsCond0Mean(:,:);
wAmpsTrigMean(:,end+1) = wAmpsBlank0Mean;

% figure
fH = [];
fH(1) = figure;
set(gcf,'Position',tsFigPos)
set(gca,'ColorOrder',trigColors)
hold all
plot(t, squeeze(wAmpsTrigMean(:,plotOrder)))
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
plot(t, mean(squeeze(wAmpsTrigMean(:,plotOrder(1:(nTrigs-1)/2))),2),'color',trigBlue,'LineWidth',4)
plot(t, mean(squeeze(wAmpsTrigMean(:,plotOrder(end-(nTrigs-1)/2):end-1)),2),'color',trigRed,'LineWidth',4)
legend(trigNames(plotOrder))
xlabel('time (ms)')
ylabel('wavelet amp')
title(sprintf('%d Hz, channel %d', ssvefFreq, channel))

% individual trials
fH(2) = figure;
set(gcf,'Position',condFigPos)
for iCue = 1:2
    for iT1 = 1:2
        for iT2 = 1:2
            iTrig = (iCue-1)*4 + (iT1-1)*2 + iT2;
            subplot(4,2,iTrig)
            plot(t, wAmpsCond0(:,:,iCue,iT1,iT2))
%             hold on
%             plot(t, 2*ones(size(t)),'k')
            title(trigNames{iTrig})
            xlim([t(1) t(end)])
            ylim([-2 4])
        end
    end
end
rd_supertitle(sprintf('%d Hz, channel %d', ssvefFreq, channel))

if saveFigs
    rd_saveAllFigs(fH, {'waveletSingleTrials','waveletSingleTrialsByCond'}, 'plot', figDir)
end

%% Hilbert on average across trials
channels = 60; % R0504_20150805: [23 43 25 13 16]; % [14 23 25 26 60]; % R0817_20150504: [13 14 23 25 43], [7 8 13 20 36];
ssvefFreq = 30;
Fbp = ssvefFreq + [-1.6 1.6];
hAmps = [];
for iTrig = 1:nTrigs
    data = trigMean(:,channels,iTrig)'; % channels by samples
    dataF = ft_preproc_bandpassfilter(data,Fs,Fbp);
    dataFH = abs(hilbert(mean(dataF,1))); % average bandpassed time series across channels
    hAmps(:,iTrig) = dataFH; 
end

fH = [];
fH(1) = figure;
set(gcf,'Position',tsFigPos)
set(gca,'ColorOrder',trigColors)
hold all
plot(t, hAmps(:,plotOrder))
legend(trigNames(plotOrder))
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
plot(t, mean(hAmps(:,plotOrder(1:(nTrigs-1)/2)),2),'color',trigBlue,'LineWidth',4)
plot(t, mean(hAmps(:,plotOrder(end-(nTrigs-1)/2):end-1),2),'color',trigRed,'LineWidth',4)
legend(trigNames(plotOrder))
xlabel('time (ms)')
ylabel('Hilbert amp')
title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels)])

% condition subplots
fH(2) = figure;
set(gcf,'Position',condFigPos)
for iTrig = 1:(nTrigs-1)/2
    subplot((nTrigs-1)/2,1,iTrig)
    set(gca,'ColorOrder',[trigBlue; trigRed])
    hold all
    plot(t, hAmps(:,iTrig*2-1:iTrig*2))
    legend(trigNames{iTrig*2-1:iTrig*2})
    for iEv = 1:numel(eventTimes)
        vline(eventTimes(iEv),'k');
    end
%     ylim([-1 2.5])
    if iTrig==1
        title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels)])
    end
end
xlabel('time (ms)')
ylabel('Hilbert amp')

% present/absent
fH(3) = figure;
set(gcf,'Position',tsFigPos)
hold on
for iTrig = 1:(nTrigs-1)/2 
    p1 = plot(t, mean(hAmps(:,iTrig*2-1:iTrig*2),2));
    set(p1, 'Color', trigColorsPA4(iTrig,:), 'LineWidth', 1.5)
end
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
legend('T1p-T2p','T1a-T2p','T1p-T2a','T1a-T2a')
xlabel('time (ms)')
ylabel('Hilbert amp')
title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels)])

% attend T1/T2 with condition error bars
mean1 = mean(hAmps(:,plotOrder(1:(nTrigs-1)/2)),2);
ste1 = std(hAmps(:,plotOrder(1:(nTrigs-1)/2)),0,2)./(sqrt((nTrigs-1)/2));
mean2 = mean(hAmps(:,plotOrder(end-(nTrigs-1)/2):end-1),2);
ste2 = std(hAmps(:,plotOrder(end-(nTrigs-1)/2):end-1),0,2)./(sqrt((nTrigs-1)/2));
fH(4) = figure;
set(gcf,'Position',tsFigPos)
hold on
shadedErrorBar(t, mean1, ste1, {'color',trigBlue,'LineWidth',4}, 1)
shadedErrorBar(t, mean2, ste2, {'color',trigRed,'LineWidth',4}, 1)
plot(t, hAmps(:,end), 'k')
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
legend('attend T1','attend T2')
xlabel('time (ms)')
ylabel('Hilbert amp')
title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels)])

if saveFigs
    figPrefix = ['plot_ch' sprintf('%d_', channels) sprintf('%dHz', ssvefFreq)];
    rd_saveAllFigs(fH, {'hilbertTrialAve','hilbertTrialAveByCond','hilbertTrialAvePA','hilbertTrialAveAttT1T2Error'}, figPrefix, figDir)
end

%% Filter single trial data
Fbp = ssvefFreq + [-1.6 1.6];
condDataF = [];
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

blankDataF = [];
for iTrial = 1:size(blankData,3)
    data = blankData(:,:,iTrial)';
    dataF = ft_preproc_bandpassfilter(data,Fs,Fbp);
    blankDataF(:,:,iTrial) = dataF';
end

% mean across trials
condDataFMean = squeeze(nanmean(condDataF,3));
blankDataFMean = squeeze(nanmean(blankDataF,3));

% let trigFMean have the conditions 1-9 in the third dimension
trigFMean = condDataFMean(:,:,:);
trigFMean(:,:,end+1) = blankDataFMean;

%% Hilbert on filtered single trial means
hTrigFMean = [];
for iChan = 1:nChannels
    for iTrig = 1:nTrigs
        dataF = trigFMean(:,iChan,iTrig);
        dataFH = abs(hilbert(dataF));
        hTrigFMean(:,iChan,iTrig) = dataFH;
    end
end

channel = 8;
figure
set(gcf,'Position',tsFigPos)
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
title(sprintf('%d Hz, channel %d', ssvefFreq, channel))

if saveFigs
    figPrefix = sprintf('plot_ch%d_%dHz', channel, ssvefFreq);
    rd_saveAllFigs(gcf, {'hilbertSingleTrials'}, figPrefix, figDir)
end

%% Time-frequency
channels = 60;
taper          = 'hanning';
foi            = 1:50;
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

tfAmpsAtt(:,:,1) = nanmean(tfAmps(:,:,plotOrder(1:(nTrigs-1)/2)),3);
tfAmpsAtt(:,:,2) = nanmean(tfAmps(:,:,plotOrder((nTrigs-1)/2+1:end-1)),3);

for iTrig = 1:(nTrigs-1)/2 
    tfAmpsPA(:,:,iTrig) = mean(tfAmps(:,:,iTrig*2-1:iTrig*2),3);
end
t1PADiff = mean(tfAmpsPA(:,:,[1 3]),3)-mean(tfAmpsPA(:,:,[2 4]),3);
t2PADiff = mean(tfAmpsPA(:,:,[1 2]),3)-mean(tfAmpsPA(:,:,[3 4]),3);

% figures
ytick = 10:10:numel(foi);
xtick = 51:50:numel(toi);
clims = [0 30];
diffClims = [-10 10];
hack = plotOrder;
hack(hack>4) = hack(hack>4)+1;

fH = [];
fH(1) = figure;
set(gcf,'Position',tf9FigPos)
for iTrig = 1:nTrigs
    subplot(2,5,hack(iTrig))
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

fH(2) = figure;
set(gcf,'Position',tf3FigPos)
attNames = {'attT1','attT2'};
for iAtt = 1:size(tfAmpsAtt,3)
    subplot(1,3,iAtt)
    imagesc(tfAmpsAtt(:,:,iAtt),clims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    xlabel('time (ms)')
    ylabel('frequency (Hz)')
    title(attNames{iAtt})
end
subplot(1,3,3)
imagesc(tfAmpsAtt(:,:,2)-tfAmpsAtt(:,:,1),diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (ms)')
ylabel('frequency (Hz)')
title('attT2 - attT1')
rd_supertitle(['channel' sprintf(' %d', channels)]);
rd_raiseAxis(gca);

fH(3) = figure;
set(gcf,'Position',tf9FigPos)
paNames = {'T1p-T2p','T1a-T2p','T1p-T2a','T1a-T2a'};
for iPA = 1:size(tfAmpsPA,3)
    subplot(2,4,iPA)
    imagesc(tfAmpsPA(:,:,iPA),clims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    xlabel('time (ms)')
    ylabel('frequency (Hz)')
    title(paNames{iPA})
end
subplot(2,4,5)
imagesc(t1PADiff,diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (ms)')
ylabel('frequency (Hz)')
title('T1 P-A')
subplot(2,4,6)
imagesc(t2PADiff,diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (ms)')
ylabel('frequency (Hz)')
title('T2 P-A')
subplot(2,4,7)
imagesc(t2PADiff - t1PADiff,diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (ms)')
ylabel('frequency (Hz)')
title('T2 vs. T1 P-A')
rd_supertitle(['channel' sprintf(' %d', channels)]);
rd_raiseAxis(gca);

if saveFigs
    figPrefix = ['im_ch' sprintf('%d_', channels)];
    rd_saveAllFigs(fH, {'timeFreqByCond','timeFreqAtt','timeFreqPA'}, figPrefix(1:end-1), figDir)
end

%% Time-frequency - single trials
channel = 60;
taper          = 'hanning';
foi            = 1:50;
t_ftimwin      = 10 ./ foi;
toi            = tstart/1000:0.01:tstop/1000;
tfSingleAmps = [];
for iTrig = 1:nTrigs-1
    [iCue,iT1,iT2] = rd_indToFactorialInd(iTrig,[2,2,2]);
    data = squeeze(condData(:,channel,:,iCue,iT1,iT2))'; % trials by samples
    [spectrum,ntaper,freqoi,timeoi] = ft_specest_mtmconvol(data, t/1000, ...
        'timeoi', toi, 'freqoi', foi, 'timwin', t_ftimwin, ...
        'taper', taper, 'dimord', 'chan_time_freqtap');
    specAmp = squeeze(nanmean(abs(spectrum),1)); % mean across trials
    tfSingleAmps(:,:,iTrig) = specAmp';
end

% blank
data = squeeze(blankData(:,channel,:))';
[spectrum,ntaper,freqoi,timeoi] = ft_specest_mtmconvol(data, t/1000, ...
    'timeoi', toi, 'freqoi', foi, 'timwin', t_ftimwin, ...
    'taper', taper, 'dimord', 'chan_time_freqtap');
specAmp = squeeze(nanmean(abs(spectrum),1)); % mean across trials
tfSingleAmps(:,:,nTrigs) = specAmp';

tfSingleAmpsAtt(:,:,1) = nanmean(tfSingleAmps(:,:,plotOrder(1:(nTrigs-1)/2)),3);
tfSingleAmpsAtt(:,:,2) = nanmean(tfSingleAmps(:,:,plotOrder((nTrigs-1)/2+1:end-1)),3);

for iTrig = 1:(nTrigs-1)/2 
    tfSingleAmpsPA(:,:,iTrig) = mean(tfSingleAmps(:,:,iTrig*2-1:iTrig*2),3);
end
t1SinglePADiff = mean(tfSingleAmpsPA(:,:,[1 3]),3)-mean(tfSingleAmpsPA(:,:,[2 4]),3);
t2SinglePADiff = mean(tfSingleAmpsPA(:,:,[1 2]),3)-mean(tfSingleAmpsPA(:,:,[3 4]),3);

% figures
ytick = 10:10:numel(foi);
xtick = 51:50:numel(toi);
clims = [0 70];
diffClims = [-10 10];
hack = plotOrder;
hack(hack>4) = hack(hack>4)+1;

fH = [];
fH(1) = figure;
set(gcf,'Position',tf9FigPos)
for iTrig = 1:nTrigs
    subplot(2,5,hack(iTrig))
    imagesc(tfSingleAmps(:,:,iTrig),clims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    if iTrig==nTrigs
        xlabel('time (ms)')
        ylabel('frequency (Hz)')
    end
    title(trigNames{iTrig})
end
rd_supertitle(sprintf('channel %d', channel));
rd_raiseAxis(gca);

fH(2) = figure;
set(gcf,'Position',tf3FigPos)
attNames = {'attT1','attT2'};
for iAtt = 1:size(tfSingleAmpsAtt,3)
    subplot(1,3,iAtt)
    imagesc(tfSingleAmpsAtt(:,:,iAtt),clims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    xlabel('time (ms)')
    ylabel('frequency (Hz)')
    title(attNames{iAtt})
end
subplot(1,3,3)
imagesc(tfSingleAmpsAtt(:,:,2)-tfSingleAmpsAtt(:,:,1),diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (ms)')
ylabel('frequency (Hz)')
title('attT2 - attT1')
rd_supertitle(sprintf('channel %d', channel));
rd_raiseAxis(gca);

fH(3) = figure;
set(gcf,'Position',tf9FigPos)
paNames = {'T1p-T2p','T1a-T2p','T1p-T2a','T1a-T2a'};
for iPA = 1:size(tfSingleAmpsPA,3)
    subplot(2,4,iPA)
    imagesc(tfSingleAmpsPA(:,:,iPA),clims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    xlabel('time (ms)')
    ylabel('frequency (Hz)')
    title(paNames{iPA})
end
subplot(2,4,5)
imagesc(t1SinglePADiff,diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (ms)')
ylabel('frequency (Hz)')
title('T1 P-A')
subplot(2,4,6)
imagesc(t2SinglePADiff,diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (ms)')
ylabel('frequency (Hz)')
title('T2 P-A')
subplot(2,4,7)
imagesc(t2SinglePADiff - t1SinglePADiff,diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (ms)')
ylabel('frequency (Hz)')
title('T2 vs. T1 P-A')
rd_supertitle(sprintf('channel %d', channel));
rd_raiseAxis(gca);

if saveFigs
    figPrefix = ['im_ch' sprintf('%d_', channel)];
    rd_saveAllFigs(fH, {'timeFreqSingleByCond','timeFreqSingleAtt','timeFreqSinglePA'}, figPrefix(1:end-1), figDir)
end

%% Bandpassed average time series
channels = 23; % [13 14 23 25 43], [7 8 13 20 36];
freq = 10;
Fbp = freq + [-1.6 1.6];
dataBP = [];
for iTrig = 1:nTrigs
    data = trigMean(:,channels,iTrig)'; % channels by samples
%     dataF = ft_preproc_bandpassfilter(data,Fs,Fbp); % does the same thing as below
%     dataBP(:,iTrig) = mean(dataF,1); 
    dataBP(:,iTrig) = ft_preproc_bandpassfilter(mean(data,1),Fs,Fbp);
end

fH = [];
fH(1) = figure;
set(gcf,'Position',tsFigPos)
set(gca,'ColorOrder',trigColors)
hold all
plot(t, dataBP(:,plotOrder))
legend(trigNames(plotOrder))
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
plot(t, mean(dataBP(:,plotOrder(1:(nTrigs-1)/2)),2),'color',trigBlue,'LineWidth',4)
plot(t, mean(dataBP(:,plotOrder(end-(nTrigs-1)/2):end-1),2),'color',trigRed,'LineWidth',4)
legend(trigNames(plotOrder))
xlabel('time (ms)')
ylabel('amplitude')
title([sprintf('%d Hz bandpass, channel', freq) sprintf(' %d', channels)])

% condition subplots
fH(2) = figure;
set(gcf,'Position',condFigPos)
for iTrig = 1:(nTrigs-1)/2
    subplot((nTrigs-1)/2,1,iTrig)
    set(gca,'ColorOrder',[trigBlue; trigRed])
    hold all
    plot(t, dataBP(:,iTrig*2-1:iTrig*2))
    legend(trigNames{iTrig*2-1:iTrig*2})
    for iEv = 1:numel(eventTimes)
        vline(eventTimes(iEv),'k');
    end
%     ylim([-1 2.5])
    if iTrig==1
        title([sprintf('%d Hz bandpass, channel', freq) sprintf(' %d', channels)])
    end
end
xlabel('time (ms)')
ylabel('amplitude')

if saveFigs
    figPrefix = ['plot_ch' sprintf('%d_', channels) sprintf('%dHz', freq)];
    rd_saveAllFigs(fH, {'bandpassTS','bandpassTSByCond'}, figPrefix, figDir)
end

%% Bandpassed single trials
condDataF = [];
for iTrial = 1:size(condData,3)
    for iCue = 1:numel(cueConds)
        for iT1 = 1:numel(t1Conds)
            for iT2 = 1:numel(t2Conds)
                data = condData(:,channels,iTrial,iCue,iT1,iT2)'; % channels by samples
                dataF = ft_preproc_bandpassfilter(mean(data,1),Fs,Fbp);
                condDataF(:,iTrial,iCue,iT1,iT2) = dataF';
            end
        end
    end
end

blankDataF = [];
for iTrial = 1:size(blankData,3)
    data = blankData(:,channels,iTrial)';
    dataF = ft_preproc_bandpassfilter(mean(data,1),Fs,Fbp);
    blankDataF(:,iTrial) = dataF';
end

fH = [];
fH(1) = figure;
set(gcf,'Position',condFigPos)
for iTrig = 1:nTrigs-1
    [iCue, iT1, iT2] = rd_indToFactorialInd(iTrig, ...
        [numel(cueConds), numel(t1Conds), numel(t2Conds)]);
    subplot(4,2,iTrig)
    p1 = plot(t, condDataF(:,:,iCue,iT1,iT2));
%     set(p1,'Color',trigColors(plotOrder(iTrig),:))
    title(trigNames{iTrig})
    xlim([t(1) t(end)])
    ylim([-500 500])
    for iEv = 1:numel(eventTimes)
        vline(eventTimes(iEv),'k');
    end
    if iTrig==7
        xlabel('time (ms)')
        ylabel('amplitude')
    end
end
rd_supertitle([sprintf('%d Hz bandpass, channel', freq) sprintf(' %d', channels)])
rd_raiseAxis(gca);

figure;
set(gcf,'Position',condFigPos)
for iTrig = 9
    subplot(4,2,1)
    p1 = plot(t, blankDataF);
%     set(p1,'Color',trigColors(plotOrder(iTrig),:))
    title(trigNames{iTrig})
    xlim([t(1) t(end)])
    ylim([-500 500])
    for iEv = 1:numel(eventTimes)
        vline(eventTimes(iEv),'k');
    end
    xlabel('time (ms)')
    ylabel('amplitude')
end


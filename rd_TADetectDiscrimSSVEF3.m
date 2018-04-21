function rd_TADetectDiscrimSSVEF3(exptDir, sessionDir, fileBase, analStr, ssvefFreq, nTopChannels, iqrThresh, weightChannels, trialSelection, respTargetSelection, exptType)

% single trial analysis

%% Setup
if nargin==0 || ~exist('exptDir','var')
    exptType = 'TANoise';
    switch exptType
        case 'TADetectDiscrim'
            exptDir = '/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG';
            sessionDir = 'R0817_20150504';
            fileBase = 'R0817_TADeDi_5.4.15';
            analStr = 'ebi'; % '', 'ebi', etc.
            ssvefFreq = 30;
            nTopChannels = 5; % 1, 5, etc., or [] for iqrThresh
            iqrThresh = []; % 10, or [] for nTopChannels
            weightChannels = 0; % weight channels according to average SSVEF amp - only works for top channels
            trialSelection = 'all'; % 'all','validCorrect', etc
            respTargetSelection = ''; % '','T1Resp','T2Resp'
            
        case 'TAContrast'
            exptDir = '/Local/Users/denison/Data/TAContrast/MEG';
            sessionDir = 'R0817_20171019';
            fileBase = 'R0817_TACont_10.19.17';
            analStr = 'ebi'; % '', 'ebi', etc.
            ssvefFreq = 20;
            nTopChannels = 5; % 1, 5, etc., or [] for iqrThresh
            iqrThresh = []; % 10, or [] for nTopChannels
            weightChannels = 0; % weight channels according to average SSVEF amp - only works for top channels
            trialSelection = 'all'; % 'all','validCorrect', etc
            respTargetSelection = ''; % '','T1Resp','T2Resp'
            
        case 'TANoise'
            exptDir = '/Local/Users/denison/Data/TANoise/MEG';
            sessionDir = 'R0817_20171212';
            fileBase = 'R0817_TANoise_12.12.17';
            analStr = 'ebi'; % '', 'ebi', etc.
            ssvefFreq = 20;
            nTopChannels = 5; % 1, 5, etc., or [] for iqrThresh
            iqrThresh = []; % 10, or [] for nTopChannels
            weightChannels = 0; % weight channels according to average SSVEF amp - only works for top channels
            trialSelection = 'all'; % 'all','validCorrect', etc
            respTargetSelection = ''; % '','T1Resp','T2Resp'
            
        otherwise
            error('exptType not recognized')
    end
end

topChannels = 1:nTopChannels;

dataDir = sprintf('%s/%s', exptDir, sessionDir);
matDir = sprintf('%s/mat', dataDir);

if ~isempty(nTopChannels) && ~isempty(iqrThresh)
    error('set either nTopChannels or iqrThresh to empty')
else
    if ~isempty(nTopChannels)
        channelSelection = 'topchannels';
        channelSelectionStr = sprintf('topChannels%d', numel(topChannels));
        if weightChannels
            channelSelectionStr = [channelSelectionStr 'W'];
        end
    elseif ~isempty(iqrThresh)
        channelSelection = 'iqrthresh';
        channelSelectionStr = sprintf('iqrThresh%d', iqrThresh);
    else
        error('set either nTopChannels or iqrThresh to a value for channel selection')
    end
end

switch analStr
    case ''
        savename = sprintf('%s/%s_ssvef_workspace.mat', matDir, fileBase);
        channelsFileName = sprintf('%s/channels_%dHz.mat', matDir, ssvefFreq);
        analysisFileName = sprintf('%s/analysis_singleTrials_%s_%s_%sTrials%s_%dHz.mat', matDir, fileBase, channelSelectionStr, trialSelection, respTargetSelection, ssvefFreq);
    otherwise
        savename = sprintf('%s/%s_%s_ssvef_workspace.mat', matDir, fileBase, analStr);
        channelsFileName = sprintf('%s/channels_%dHz_%s.mat', matDir, ssvefFreq, analStr);
        analysisFileName = sprintf('%s/analysis_singleTrials_%s_%s_%s_%sTrials%s_%dHz.mat', matDir, fileBase, analStr, channelSelectionStr, trialSelection, respTargetSelection, ssvefFreq);
end

%% Get the data
load(savename)

%% Update behav
behav = behavior(behav);

%% Settings after loading the data
saveAnalysis = 0;
saveFigs = 0;
plotFigs = 1;

excludeTrialsFt = 1;
excludeSaturatedEpochs = 1;

load(channelsFileName);
switch channelSelection
    case 'topchannels'
        channels = channelsRanked(topChannels);
    case 'iqrthresh'
        channels = find(peakMeansStimAve > median(peakMeansBlank) + iqrThresh*iqr(peakMeansBlank));
    otherwise
        error('channelSelection not recognized')
end
if isempty(channels)
    fprintf('No channels found for %s iqrThresh %d Hz ... exiting.', sessionDir, ssvefFreq)
    return
end

% channel weights
if weightChannels && strcmp(channelSelection,'topchannels')
    chw = (channelsRankedAmps(topChannels)/channelsRankedAmps(1))';
    wstr = 'W';
    wstr2 = 'W_';
    wstrt = ' W';
else
    chw = ones(size(topChannels))';
    wstr = '';
    wstr2 = '';
    wstrt = '';
end

%% Store settings for this analysis
A.fileBase = fileBase;
A.analStr = analStr;
A.excludeTrialsFt = excludeTrialsFt;
A.excludeSaturatedEpochs = excludeSaturatedEpochs;
A.ssvefFreq = ssvefFreq;
A.channels = channels;
A.chw = chw;
A.Fs = Fs;
A.t = t;
if ~exist('eventTimes','var')
    eventTimes = [0 500 1500 2100 3100];
end
A.eventTimes = eventTimes;
A.trigNames = trigNames;

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
    trigData(:,:,trials_rejected) = NaN;
    
    % update figDir
    figDir = [figDir '_ft'];
    
    % update analysis file
    switch analStr
        case ''
            analysisFileName = sprintf('%s/analysis_singleTrials_%s_ft_%s_%sTrials%s_%dHz.mat', matDir, fileBase, channelSelectionStr, trialSelection, respTargetSelection, ssvefFreq);
        otherwise
            analysisFileName = sprintf('%s/analysis_singleTrials_%s_%s_ft_%s_%sTrials%s_%dHz.mat', matDir, fileBase, analStr, channelSelectionStr, trialSelection, respTargetSelection, ssvefFreq);
    end
end

%% Make figDir if needed
figDir = sprintf('%s_singleTrials_%s_%sTrials%s', figDir, channelSelectionStr, trialSelection, respTargetSelection);

if ~exist(figDir,'dir') && saveFigs
    mkdir(figDir)
end

%% Organize trials into conditions
switch exptType
    case 'TADetectDiscrim'
        targetCondNames = {'target type T1','target type T2'};
        t1Conds = {[1 2], 0}; % present, absent
        t2Conds = {[1 2], 0}; % present, absent
    case {'TAContrast','TANoise'}
        targetCondNames = {'target pedestal T1','target pedestal T2'};
        t1Conds = {1, 2}; % pedestal decrement, pedestal increment
        t2Conds = {1, 2}; % pedestal decrement, pedestal increment
    otherwise
        error('exptType not recognized')
end

cueCondIdx = strcmp(behav.responseData_labels, 'cue condition');
t1CondIdx = strcmp(behav.responseData_labels, targetCondNames{1});
t2CondIdx = strcmp(behav.responseData_labels, targetCondNames{2});
nTrials = size(behav.responseData_all,1);

blankCond = 1;
cueConds = {[2 3], [4 5]}; % cue T1, cue T2

switch respTargetSelection
    case 'T1Resp'
        rSelect = behav.responseTarget==1;
    case 'T2Resp'
        rSelect = behav.responseTarget==2;
    case ''
        rSelect = ones(nTrials,1);
    otherwise
        error('respTargetSelection not recognized')
end

switch trialSelection
    case 'correct'
        wSelect = behav.acc==1;
    case 'incorrect'
        wSelect = behav.acc==0;
    case 'validCorrect'
        wValid = behav.cueValidity==1;
        wCorrect = behav.acc==1; 
        wSelect = wValid & wCorrect;
    case 'detectHit'
        wSelect = behav.detectHMFC(:,1)==1;
    case 'detectMiss'
        wSelect = behav.detectHMFC(:,2)==1;
    case 'detectFA'
        wSelect = behav.detectHMFC(:,3)==1;
    case 'detectCR'
        wSelect = behav.detectHMFC(:,4)==1;
    case 'discrimCorrect'
        wSelect = behav.discrimCI(:,1)==1;
    case 'discrimIncorrect'
        wSelect = behav.discrimCI(:,2)==1;
    case 'all'
        wSelect = ones(nTrials,1);
    otherwise
        error('trialSelection not recognized')
end
wSelect = wSelect & rSelect;

trigDataSelected = trigData; % make a copy so we use it for condData but not blankData
trigDataSelected(:,:,wSelect~=1)=NaN;

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
            nTrialsCond(iCue,iT1,iT2) = nnz(w & wSelect);
            fprintf('Number of trials %d %d %d: %d\n', iCue, iT1, iT2, nnz(w & wSelect))
            
            condData(:,:,:,iCue,iT1,iT2) = trigDataSelected(:,:,w);
        end
    end
end

wBlank = behav.responseData_all(:,cueCondIdx) == blankCond;
blankData = trigData(:,:,wBlank);

% let trigMean have the conditions 1-9 in the fourth dimension
%%% note that here we are selecting channels!
%%% and that trigMean is not a mean! we're just keeping the name.
trigMean = reshape(condData(:,channels,:,:,:,:), ...
    size(condData,1), numel(channels), size(condData,3), ...
    numel(cueConds)*numel(t1Conds)*numel(t2Conds));

% randomly sample trials from blank to match number in other conditions
randomBlankTrials = randperm(size(blankData,3));
trigMean(:,:,:,end+1) = blankData(:,channels,randomBlankTrials(1:size(condData,3)));
nTrigs = size(trigMean,4);
nTrialsPerCond = size(trigMean,3);

A.nTrialsCond = nTrialsCond;
A.trigMean = trigMean;

%% FFT on mean time series for each trigger type
% do the fft for each channel
% nfft = 2^nextpow2(nSamples); % Next power of 2 from length of y
% Y = fft(trigMean,nfft)/nSamples; % Scale by number of samples

% only go from cue to post-cue
tidx1 = find(t==eventTimes(2));
tidx2 = find(t==eventTimes(5))-1;
nfft = numel(tidx1:tidx2);
Y = fft(trigMean(tidx1:tidx2,:,:,:),nfft)/nfft; % Scale by number of samples
f = Fs/2*linspace(0,1,nfft/2+1); % Fs/2 is the maximum frequency that can be measured
amps = 2*abs(Y(1:nfft/2+1,:,:,:)); % Multiply by 2 since only half the energy is in the positive half of the spectrum?

A.f = f;
A.Y = Y;
A.amps = amps;

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
tf3FigPos = [200 475 980 330];
paau3FigPos = [800 30 450 830];
paau6FigPos = [0 90 1000 800];
tf9SquareFigPos = [50 50 850 850];
tf6SquareFigPos = [50 50 850 530];

set(0,'defaultLineLineWidth',1)

%% Time series and FFT
% mean across channels
trigMeanMean = squeeze(rd_wmean(trigMean,chw,2));
ampsMean = squeeze(rd_wmean(amps,chw,2));
A.trigMeanMean = trigMeanMean;
A.ampsMean = ampsMean;

if plotFigs
fH = [];
fH(1) = figure;
set(gcf,'Position',ts2FigPos)

% time
subplot(3,1,1)
hold on
for iTrig = 1:nTrigs
    plot(t, squeeze(nanmean(trigMeanMean(:,:,plotOrder(iTrig)),2)),'color',trigColors(iTrig,:))
end
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
xlim([t(1) t(end)])
xlabel('time (ms)')
ylabel('amplitude')
title(['channel' sprintf(' %d', channels) wstrt])

% frequency
subplot(3,1,2)
hold on
for iTrig = 1:nTrigs
    plot(f, squeeze(nanmean(ampsMean(:,:,plotOrder(iTrig)),2)),'color',trigColors(iTrig,:))
end
xlim([1 200])
ylim([0 60])
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
legend(trigNames(plotOrder))

subplot(3,1,3)
hold on
plot(f, squeeze(nanmean(mean(ampsMean(:,:,1:end-1),3),2)),'color',[.66 .5 .78])
plot(f, squeeze(nanmean(ampsMean(:,:,end),2)),'color',trigColors(end,:))
xlim([1 200])
ylim([0 60])
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
legend('stim average','blank')

% log scale
fH(2) = figure;
set(gcf,'Position',tsFigPos)
loglog(f, squeeze(nanmean(mean(ampsMean(:,:,1:end-1),3),2)),'color',[.66 .5 .78])
hold on
loglog(f, squeeze(nanmean(ampsMean(:,:,end),2)),'color',trigColors(end,:))
vline(ssvefFreq,'color','k','LineStyle',':')
xlim([0 200])
xlabel('Log frequency (Hz)')
ylabel('log(|Y(f)|)')
legend('stim average','blank')
legend boxoff
title(['channel' sprintf(' %d', channels) wstrt])
end

if saveFigs
    if numel(channels)==1
        figPrefix = sprintf('plot_ch%d', channels);
    else
        figPrefix = ['plot_ch' sprintf('%d_', channels(1:end-1)) sprintf('%d', channels(end)) wstr];
    end
    rd_saveAllFigs(fH, {'tsFFT','FFTLog'}, figPrefix, figDir)
end

%% Target present vs. absent
targetWindow = [-100 600];
t1Tidx = t>=eventTimes(3) + targetWindow(1) & t<=eventTimes(3) + targetWindow(2);
t2Tidx = t>=eventTimes(4) + targetWindow(1) & t<=eventTimes(4) + targetWindow(2);

% calculate pres/abs x att/unattend for each target
paauT(:,:,1,1) = [trigMeanMean(t1Tidx,:,1) trigMeanMean(t1Tidx,:,5)]; % present/attended
paauT(:,:,2,1) = [trigMeanMean(t1Tidx,:,2) trigMeanMean(t1Tidx,:,6)]; % present/unattended
paauT(:,:,3,1) = [trigMeanMean(t1Tidx,:,3) trigMeanMean(t1Tidx,:,7)] ; % absent/attended
paauT(:,:,4,1) = [trigMeanMean(t1Tidx,:,4) trigMeanMean(t1Tidx,:,8)]; % absent/unattended

paauT(:,:,1,2) = [trigMeanMean(t2Tidx,:,2) trigMeanMean(t2Tidx,:,4)];
paauT(:,:,2,2) = [trigMeanMean(t2Tidx,:,1) trigMeanMean(t2Tidx,:,3)];
paauT(:,:,3,2) = [trigMeanMean(t2Tidx,:,6) trigMeanMean(t2Tidx,:,8)];
paauT(:,:,4,2) = [trigMeanMean(t2Tidx,:,5) trigMeanMean(t2Tidx,:,7)];

% ave pres - ave abs. trial pairing not meaningful
targetPADiff = squeeze((paauT(:,:,1,:)+paauT(:,:,2,:))/2 - ...
    (paauT(:,:,3,:)+paauT(:,:,4,:))/2); 
targetPADiffAll = targetPADiff(:,:);

% calculate power only after the time of the target presentation
tZeroIdx = find(targetWindow(1):targetWindow(2)==0);
targetNfft = 2^nextpow2(targetWindow(2)+1);
targetY = fft(targetPADiff(tZeroIdx:end,:,:),targetNfft)/(targetWindow(2)+1);
targetF = Fs/2*linspace(0,1,targetNfft/2+1);
targetAmps = 2*abs(targetY(1:targetNfft/2+1,:,:));
targetAmpsAll = targetAmps(:,:);

% store results
A.targetWindow = targetWindow;
A.paauT = paauT;
A.targetPADiff = targetPADiff;
A.targetF = targetF;
A.targetPADiffAmps = targetAmps;

if plotFigs
switch exptType
    case 'TADetectDiscrim'
        names = {'target present','target absent'};
    case 'TAContrast'
        names = {'target decrement','target increment'};
    case 'TANoise'
        names = {'vertical','horizontal'};
end
colors = get(gca,'ColorOrder');
fH = [];
fH(1) = figure;
set(gcf,'Position',ts3FigPos)
subplot(3,1,1)
hold on
plot(targetWindow(1):targetWindow(2), squeeze(nanmean(targetPADiff,2)))
legend('T1','T2')
[~, emp, err] = rd_bootstrapCI(targetPADiff(:,:,1)');
shadedErrorBar(targetWindow(1):targetWindow(2), emp, err, {'color',colors(1,:),'LineWidth',2}, 1)
[~, emp, err] = rd_bootstrapCI(targetPADiff(:,:,2)');
shadedErrorBar(targetWindow(1):targetWindow(2), emp, err, {'color',colors(2,:),'LineWidth',2}, 1)
vline(0,'k')
xlabel('time (ms)')
ylabel('\Delta amplitude')
title(sprintf('%s - %s', names{1}, names{2}))
subplot(3,1,2)
hold on
[~, emp, err] = rd_bootstrapCI(targetPADiffAll');
shadedErrorBar(targetWindow(1):targetWindow(2), emp, err, {'color','k','LineWidth',2}, 1)
vline(0,'k')
legend('T1 & T2')
xlabel('time (ms)')
ylabel('\Delta amplitude')
subplot(3,1,3)
hold on
[~, emp, err] = rd_bootstrapCI(targetAmpsAll');
shadedErrorBar(targetF, emp, err, {'color','k','LineWidth',2}, 1)
plot(targetF, nanmean(targetAmpsAll,2), 'b')
xlim([0 50])
xlabel('frequency (Hz)')
ylabel('amplitude')
rd_supertitle2(['channel' sprintf(' %d', channels) wstrt])
end

if saveFigs
    if numel(channels)==1
        figPrefix = sprintf('plot_ch%d', channels);
    else
        figPrefix = ['plot_ch' sprintf('%d_', channels(1:end-1)) sprintf('%d', channels(end)) wstr];
    end
    rd_saveAllFigs(fH, {'targetPATrialAveDiff'}, figPrefix, figDir)
end

%% Wavelet
switch ssvefFreq
    case 11
        width = 4;
    case 15
        width = 6;
    case {20, 25}
        width = 8;
    case 30
        width = 12; % 12 for 30 Hz, 16 for 40 Hz gives 127 ms duration, 5 Hz bandwidth
    case 40
        width = 16;
    otherwise
        error('ssvefFreq not recognized')
end
wBaselineWindow = NaN;
% wBaselineWindow = [-500 0]; % [-300 -200];
% wBaselineWindowIdx = find(t==wBaselineWindow(1)):find(t==wBaselineWindow(2));

% only frequency of interest
wAmps0 = [];
wITPC0 = [];
wSpec0 = [];
foi = ssvefFreq;
for iTrig = 1:nTrigs
    for iCh = 1:numel(channels)
        data = squeeze(trigMean(:,iCh,:,iTrig))'; % trials by samples
        [spectrum,freqoi,timeoi] = ft_specest_wavelet(data, t/1000, 'freqoi', foi, 'width', width);
        spec = squeeze(spectrum);
        specAmp = abs(squeeze(spectrum));
        itpc = squeeze(abs(nanmean(exp(1i*angle(spectrum)),1))); % mean across trials
        
        if all(size(specAmp)>1) % if two-dimensional
            wAmp = specAmp;
            spec = spec;
        else
            wAmp = specAmp';
            spec = spec';
        end
        %     wAmpNorm = wAmp./nanmean(nanmean(wAmp(:,wBaselineWindowIdx)))-1;
        %     wAmps0(:,:,iTrig) = wAmpNorm';
        wAmps0(:,iCh,:,iTrig) = wAmp';
        wITPC0(:,iCh,iTrig) = itpc;
        wSpec0(:,:,iTrig,iCh) = spec';
    end
end
wAmps = squeeze(rd_wmean(wAmps0,chw,2)); % mean across channels
wITPC = squeeze(rd_wmean(wITPC0,chw,2));

% attT1T2 combined
att1 = []; att2 = [];
att1Spec = []; att2Spec = [];
conds1 = plotOrder(1:(nTrigs-1)/2);
conds2 = plotOrder((nTrigs-1)/2+1:nTrigs-1);
for i=1:4
    att1 = cat(2, att1, wAmps(:,:,conds1(i)));
    att2 = cat(2, att2, wAmps(:,:,conds2(i)));
    att1Spec = cat(2, att1Spec, wSpec0(:,:,conds1(i),:)); 
    att2Spec = cat(2, att2Spec, wSpec0(:,:,conds2(i),:)); 
end
wAmpsAtt(:,:,1) = att1;
wAmpsAtt(:,:,2) = att2;
wSpecAtt(:,:,1,:) = att1Spec; % time x trials x att cond x channels
wSpecAtt(:,:,2,:) = att2Spec;
attNames = {'attT1','attT2'};

% PA combined
wAmpsPA(:,:,1) = [wAmps(:,:,1) wAmps(:,:,2)];
wAmpsPA(:,:,2) = [wAmps(:,:,3) wAmps(:,:,4)];
wAmpsPA(:,:,3) = [wAmps(:,:,5) wAmps(:,:,6)];
wAmpsPA(:,:,4) = [wAmps(:,:,7) wAmps(:,:,8)];

wSpecPA(:,:,1,:) = cat(2, wSpec0(:,:,1,:), wSpec0(:,:,2,:));
wSpecPA(:,:,2,:) = cat(2, wSpec0(:,:,3,:), wSpec0(:,:,4,:));
wSpecPA(:,:,3,:) = cat(2, wSpec0(:,:,5,:), wSpec0(:,:,6,:));
wSpecPA(:,:,4,:) = cat(2, wSpec0(:,:,7,:), wSpec0(:,:,8,:));

switch exptType
    case 'TADetectDiscrim'
        PANames = {'T1p-T2p','T1a-T2p','T1p-T2a','T1a-T2a'};
        PADiffNames = 'P-A';
        xtickint = 50;
    case 'TAContrast'
        PANames = {'T1d-T2d','T1i-T2d','T1d-T2i','T1i-T2i'};
        PADiffNames = 'D-I';
        xtickint = 100;
    case 'TANoise'
        PANames = {'T1v-T2v','T1h-T2v','T1v-T2h','T1h-T2h'};
        PADiffNames = 'V-H';
        xtickint = 100;
end

% All combined
wAmpsAll = [wAmpsAtt(:,:,1) wAmpsAtt(:,:,2)];
wSpecAll = squeeze(cat(2, wSpecAtt(:,:,1,:), wSpecAtt(:,:,2,:)));

% ITPC
% Method 1: we're taking the means across conditions rather than
% recomputing the itpc for all trials in a given attention condition. can
% reconsider this choice.
% attT1T2 combined
% conds1 = plotOrder(1:(nTrigs-1)/2);
% conds2 = plotOrder((nTrigs-1)/2+1:nTrigs-1);
% wITPCAtt(:,1) = nanmean(wITPC(:,conds1),2);
% wITPCAtt(:,2) = nanmean(wITPC(:,conds2),2);
% 
% % PA combined
% wITPCPA(:,1) = nanmean(wITPC(:,[1 2]),2);
% wITPCPA(:,2) = nanmean(wITPC(:,[3 4]),2);
% wITPCPA(:,3) = nanmean(wITPC(:,[5 6]),2);
% wITPCPA(:,4) = nanmean(wITPC(:,[7 8]),2); 
% 
% % All combined
% wITPCAll = nanmean(wITPC,2);

% Method 2: compute itpc from all trials in condition
% attT1T2 combined
for iAtt = 1:numel(attNames)
    for iCh = 1:numel(channels)
        spectrum = wSpecAtt(:,:,iAtt,iCh)';
        itpc = squeeze(abs(nanmean(exp(1i*angle(spectrum)),1))); % mean across trials
        
        wITPCAtt0(:,iCh,iAtt) = itpc;
    end
end
wITPCAtt = squeeze(rd_wmean(wITPCAtt0,chw,2)); % mean across channels

% PA combined
for iPA = 1:numel(PANames)
    for iCh = 1:numel(channels)
        spectrum = wSpecPA(:,:,iPA,iCh)';
        itpc = squeeze(abs(nanmean(exp(1i*angle(spectrum)),1))); % mean across trials
        
        wITPCPA0(:,iCh,iPA) = itpc;
    end
end
wITPCPA = squeeze(rd_wmean(wITPCPA0,chw,2));

% All combined
for iCh = 1:numel(channels)
    spectrum = wSpecAll(:,:,iCh)';
    itpc = squeeze(abs(nanmean(exp(1i*angle(spectrum)),1))); % mean across trials
    
    wITPCAll0(:,iCh) = itpc;
end
wITPCAll = squeeze(rd_wmean(wITPCAll0,chw,2));


% store results
A.attNames = attNames;
A.PANames = PANames;
A.PADiffNames = PADiffNames;
A.wBaselineWindow = wBaselineWindow;
A.wAmps = wAmps;
A.wAmpsAtt = wAmpsAtt;
A.wAmpsPA = wAmpsPA;
A.wAmpsAll = wAmpsAll;

A.wSpec = wSpec0;
A.wSpecAtt = wSpecAtt;
A.wSpecPA = wSpecPA;
A.wSpecAll = wSpecAll;

A.wITPC = wITPC;
A.wITPCAtt = wITPCAtt;
A.wITPCPA = wITPCPA;
A.wITPCAll = wITPCAll;

if plotFigs
fH = [];
fH(1) = figure;
set(gcf,'Position',tsFigPos)
set(gca,'ColorOrder',trigColors)
hold all
plot(t, squeeze(nanmean(wAmps(:,:,plotOrder),2)))
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
legend(trigNames(plotOrder))
xlabel('time (ms)')
ylabel('wavelet amp')
title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels) wstrt])

fH(2) = figure;
set(gcf,'Position',tsFigPos)
hold on
plot(t, nanmean(wAmpsAtt(:,:,1),2),'color',trigBlue,'LineWidth',4)
plot(t, nanmean(wAmpsAtt(:,:,2),2),'color',trigRed,'LineWidth',4)
legend(attNames)
[~, emp, err] = rd_bootstrapCI(wAmpsAtt(:,:,1)');
shadedErrorBar(t, emp, err, {'color',trigBlue,'LineWidth',4}, 1)
[~, emp, err] = rd_bootstrapCI(wAmpsAtt(:,:,2)');
shadedErrorBar(t, emp, err, {'color',trigRed,'LineWidth',4}, 1)
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
xlabel('time (ms)')
ylabel('wavelet amp')
title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels) wstrt])

% condition subplots
fH(3) = figure;
set(gcf,'Position',condFigPos)
for iTrig = 1:(nTrigs-1)/2
    subplot((nTrigs-1)/2,1,iTrig)
    hold on
    plot(t, nanmean(wAmps(:,:,iTrig*2-1),2),'color',trigBlue,'LineWidth',4)
    plot(t, nanmean(wAmps(:,:,iTrig*2),2),'color',trigRed,'LineWidth',4)
    legend(trigNames{iTrig*2-1:iTrig*2})
    [~, emp, err] = rd_bootstrapCI(wAmps(:,:,iTrig*2-1)');
    shadedErrorBar(t, emp, err, {'color',trigBlue,'LineWidth',4}, 1)
    [~, emp, err] = rd_bootstrapCI(wAmps(:,:,iTrig*2)');
    shadedErrorBar(t, emp, err, {'color',trigRed,'LineWidth',4}, 1)
%     ylim([-1 2.5])
    for iEv = 1:numel(eventTimes)
        vline(eventTimes(iEv),'k');
    end
    if iTrig==1
        title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels) wstrt])
    end
end
xlabel('time (ms)')
ylabel('wavelet amp')

% present/absent
fH(4) = figure;
set(gcf,'Position',tsFigPos)
hold on
for iPA = 1:4
    p1 = plot(t, nanmean(wAmpsPA(:,:,iPA),2));
    set(p1, 'Color', trigColorsPA4(iPA,:), 'LineWidth', 4)
end
legend(PANames)
for iPA = 1:4
    [~, emp, err] = rd_bootstrapCI(wAmpsPA(:,:,iPA)');
    shadedErrorBar(t, emp, err, {'color',trigColorsPA4(iPA,:),'LineWidth',4}, 1) 
end
% ylim([-1 2.5])
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
xlabel('time (ms)')
ylabel('wavelet amp')
title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels) wstrt])

fH(5) = figure;
set(gcf,'Position',tsFigPos)
[~, emp, err] = rd_bootstrapCI(wAmpsAll');
shadedErrorBar(t, emp, err, {'color','k','LineWidth',4})
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
legend('all trials')
xlabel('time (ms)')
ylabel('wavelet amp')
title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels) wstrt])
end

if saveFigs
    figPrefix = ['plot_ch' sprintf('%d_', channels) wstr2 sprintf('%dHz', ssvefFreq)];
    rd_saveAllFigs(fH, {'waveletTrialAve','waveletTrialAveAtt','waveletTrialAveByCond','waveletTrialAvePA','waveletTrialAveAll'}, figPrefix, figDir)
end

%% PAAU
twin = [-600 600];
t1Tidx = find(t==eventTimes(3)+twin(1)):find(t==eventTimes(3)+twin(2));
t2Tidx = find(t==eventTimes(4)+twin(1)):find(t==eventTimes(4)+twin(2));

% calculate pres/abs x att/unattend for each target
wPAAUT(:,:,1,1) = [wAmps(t1Tidx,:,1) wAmps(t1Tidx,:,5)]; % present/attended
wPAAUT(:,:,2,1) = [wAmps(t1Tidx,:,2) wAmps(t1Tidx,:,6)]; % present/unattended
wPAAUT(:,:,3,1) = [wAmps(t1Tidx,:,3) wAmps(t1Tidx,:,7)] ; % absent/attended
wPAAUT(:,:,4,1) = [wAmps(t1Tidx,:,4) wAmps(t1Tidx,:,8)]; % absent/unattended

wPAAUT(:,:,1,2) = [wAmps(t2Tidx,:,2) wAmps(t2Tidx,:,4)];
wPAAUT(:,:,2,2) = [wAmps(t2Tidx,:,1) wAmps(t2Tidx,:,3)];
wPAAUT(:,:,3,2) = [wAmps(t2Tidx,:,6) wAmps(t2Tidx,:,8)];
wPAAUT(:,:,4,2) = [wAmps(t2Tidx,:,5) wAmps(t2Tidx,:,7)];

% present vs. absent and attended vs. unattended
for iT = 1:2
    wPAT(:,:,1,iT) = cat(2, wPAAUT(:,:,1,iT), wPAAUT(:,:,2,iT)); % present
    wPAT(:,:,2,iT) = cat(2, wPAAUT(:,:,3,iT), wPAAUT(:,:,4,iT)); % absent

    wAUT(:,:,1,iT) = cat(2, wPAAUT(:,:,1,iT), wPAAUT(:,:,3,iT)); % attended
    wAUT(:,:,2,iT) = cat(2, wPAAUT(:,:,2,iT), wPAAUT(:,:,4,iT)); % unattended 
end

% combining across T1 and T2
for iPAAU = 1:4
    wPAAU(:,:,iPAAU) = cat(2, wPAAUT(:,:,iPAAU,1), wPAAUT(:,:,iPAAU,2));
end

wPA(:,:,1) = cat(2, wPAAU(:,:,1), wPAAU(:,:,2)); % present
wPA(:,:,2) = cat(2, wPAAU(:,:,3), wPAAU(:,:,4)); % absent

wAU(:,:,1) = cat(2, wPAAU(:,:,1), wPAAU(:,:,3)); % attended
wAU(:,:,2) = cat(2, wPAAU(:,:,2), wPAAU(:,:,4)); % unattended

% store results
A.wtwin = twin;
A.wt1Tidx = t1Tidx;
A.wt2Tidx = t2Tidx;
A.wPAAUT = wPAAUT;
A.wPAT = wPAT;
A.wAUT = wAUT;
A.wPAAU = wPAAU;
A.wPA = wPA;
A.wAU = wAU;

switch exptType
    case 'TADetectDiscrim'
        PAAUNames = {'P-att','P-unatt','A-att','A-unatt'};
    case 'TAContrast'
        PAAUNames = {'D-att','D-unatt','I-att','I-unatt'};
    case 'TANoise'
        PAAUNames = {'V-att','V-unatt','H-att','H-unatt'};
end

if plotFigs
% separate T1 and T2
fH = [];
fH(1) = figure;
set(gcf,'Position',paau6FigPos)
for iT = 1:2
    colors = get(gca,'ColorOrder');
    subplot(3,2,iT)
    hold on
    for iPAAU = 1:4
        p1 = plot(twin(1):twin(end), nanmean(wPAAUT(:,:,iPAAU,iT),2));
        set(p1, 'Color', colors(iPAAU,:), 'LineWidth', 4)
    end
    if iT==2
        legend(PAAUNames)
    end
    for iPAAU = 1:4
        [~, emp, err] = rd_bootstrapCI(wPAAUT(:,:,iPAAU,iT)');
        shadedErrorBar(twin(1):twin(end), emp, err, {'color',colors(iPAAU,:),'LineWidth',4}, 1)
    end
    vline(0,'k');
    xlabel('time (ms)')
    ylabel('wavelet amp')
    title(sprintf('T%d',iT))
    
    colors = [trigBlue; trigRed];
    subplot(3,2,iT+2)
    hold on
    for iAU = 1:2
        p1 = plot(twin(1):twin(end), nanmean(wAUT(:,:,iAU,iT),2));
        set(p1, 'Color', colors(iAU,:), 'LineWidth', 4)
    end
    if iT==2
        legend('att','unatt')
    end
    for iAU = 1:2
        [~, emp, err] = rd_bootstrapCI(wAUT(:,:,iAU,iT)');
        shadedErrorBar(twin(1):twin(end), emp, err, {'color',colors(iAU,:),'LineWidth',4}, 1)
    end
    vline(0,'k');
    xlabel('time (ms)')
    ylabel('wavelet amp')
    title(sprintf('T%d',iT))
    
    colors = trigColorsPA4([1 4],:);
    subplot(3,2,iT+4)
    hold on
    for iPA = 1:2
        p1 = plot(twin(1):twin(end), nanmean(wPAT(:,:,iPA,iT),2));
        set(p1, 'Color', colors(iPA,:), 'LineWidth', 4)
    end
    if iT==2
        legend(names)
    end
    for iPA = 1:2
        [~, emp, err] = rd_bootstrapCI(wPAT(:,:,iPA,iT)');
        shadedErrorBar(twin(1):twin(end), emp, err, {'color',colors(iPA,:),'LineWidth',4}, 1)
    end
    vline(0,'k');
    xlabel('time (ms)')
    ylabel('wavelet amp')
    title(sprintf('T%d',iT))
end
rd_supertitle2([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels) wstrt])

% combined across T1 and T2
fH(2) = figure;
set(gcf,'Position',paau3FigPos)
colors = get(gca,'ColorOrder');
subplot(3,1,1)
hold on
for iPAAU = 1:4
    p1 = plot(twin(1):twin(end), nanmean(wPAAU(:,:,iPAAU),2));
    set(p1, 'Color', colors(iPAAU,:), 'LineWidth', 4)
end
legend(PAAUNames)
for iPAAU = 1:4
    [~, emp, err] = rd_bootstrapCI(wPAAU(:,:,iPAAU)');
    shadedErrorBar(twin(1):twin(end), emp, err, {'color',colors(iPAAU,:),'LineWidth',4}, 1)
end
vline(0,'k');
xlabel('time (ms)')
ylabel('wavelet amp')
title('T1 & T2')

colors = [trigBlue; trigRed];
subplot(3,1,2)
hold on
for iAU = 1:2
    p1 = plot(twin(1):twin(end), nanmean(wAU(:,:,iAU),2));
    set(p1, 'Color', colors(iAU,:), 'LineWidth', 4)
end
legend('att','unatt')
for iAU = 1:2
    [~, emp, err] = rd_bootstrapCI(wAU(:,:,iAU)');
    shadedErrorBar(twin(1):twin(end), emp, err, {'color',colors(iAU,:),'LineWidth',4}, 1)
end
vline(0,'k');
xlabel('time (ms)')
ylabel('wavelet amp')
title('T1 & T2')

colors = trigColorsPA4([1 4],:);
subplot(3,1,3)
hold on
for iPA = 1:2
    p1 = plot(twin(1):twin(end), nanmean(wPA(:,:,iPA),2));
    set(p1, 'Color', colors(iPA,:), 'LineWidth', 4)
end
legend(names)
for iPA = 1:2
    [~, emp, err] = rd_bootstrapCI(wPA(:,:,iPA)');
    shadedErrorBar(twin(1):twin(end), emp, err, {'color',colors(iPA,:),'LineWidth',4}, 1)
end
vline(0,'k');
xlabel('time (ms)')
ylabel('wavelet amp')
title('T1 & T2')

rd_supertitle2([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels) wstrt])
end

if saveFigs
    figPrefix = ['plot_ch' sprintf('%d_', channels) wstr2 sprintf('%dHz', ssvefFreq)];
    rd_saveAllFigs(fH, {'waveletPAAU','waveletPAAUT1T2Comb'}, figPrefix, figDir)
end

%% ITPC plots
if plotFigs
fH = [];
fH(1) = figure;
set(gcf,'Position',tsFigPos)
set(gca,'ColorOrder',trigColors)
hold all
plot(t, wITPC(:,plotOrder))
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
legend(trigNames(plotOrder))
xlabel('time (ms)')
ylabel('wavelet itpc')
title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels) wstrt])

fH(2) = figure;
set(gcf,'Position',tsFigPos)
hold on
plot(t, wITPCAtt(:,1),'color',trigBlue,'LineWidth',2)
plot(t, wITPCAtt(:,2),'color',trigRed,'LineWidth',2)
legend(attNames)
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
xlabel('time (ms)')
ylabel('wavelet itpc')
title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels) wstrt])

% condition subplots
fH(3) = figure;
set(gcf,'Position',condFigPos)
for iTrig = 1:(nTrigs-1)/2
    subplot((nTrigs-1)/2,1,iTrig)
    hold on
    plot(t, wITPC(:,iTrig*2-1),'color',trigBlue,'LineWidth',2)
    plot(t, wITPC(:,iTrig*2),'color',trigRed,'LineWidth',2)
    legend(trigNames{iTrig*2-1:iTrig*2})
%     ylim([-1 2.5])
    for iEv = 1:numel(eventTimes)
        vline(eventTimes(iEv),'k');
    end
    if iTrig==1
        title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels) wstrt])
    end
end
xlabel('time (ms)')
ylabel('wavelet itpc')

% present/absent
fH(4) = figure;
set(gcf,'Position',tsFigPos)
hold on
for iPA = 1:4
    p1 = plot(t, wITPCPA(:,iPA));
    set(p1, 'Color', trigColorsPA4(iPA,:), 'LineWidth', 2)
end
legend(PANames)
% ylim([-1 2.5])
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
xlabel('time (ms)')
ylabel('wavelet itpc')
title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels) wstrt])

% present/absent
fH(5) = figure;
set(gcf,'Position',tsFigPos)
p1 = plot(t, wITPCAll);
set(p1, 'Color', 'k', 'LineWidth', 2)
legend('all trials')
% ylim([-1 2.5])
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
xlabel('time (ms)')
ylabel('wavelet itpc')
title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels) wstrt])
end

if saveFigs
    figPrefix = ['plot_ch' sprintf('%d_', channels) wstr2 sprintf('%dHz', ssvefFreq)];
    rd_saveAllFigs(fH, {'itpc','itpcAtt','itpcByCond','itpcPA','itpcAll'}, figPrefix, figDir)
end


%% Time-frequency - single trials
taper          = 'hanning';
foi            = 1:50;
t_ftimwin      = 10 ./ foi;
toi            = tstart/1000:0.01:tstop/1000;
tfSingleAmps0 = [];
tfSingleITPC0 = [];
for iCh = 1:numel(channels)
    channel = channels(iCh);
    for iTrig = 1:nTrigs-1
        [iCue,iT1,iT2] = rd_indToFactorialInd(iTrig,[2,2,2]);
        data = squeeze(condData(:,channel,:,iCue,iT1,iT2))'; % trials by samples
        [spectrum,ntaper,freqoi,timeoi] = ft_specest_mtmconvol(data, t/1000, ...
            'timeoi', toi, 'freqoi', foi, 'timwin', t_ftimwin, ...
            'taper', taper, 'dimord', 'chan_time_freqtap');
        specAmp = squeeze(nanmean(abs(spectrum),1)); % mean across trials
        itpc = squeeze(abs(nanmean(exp(1i*angle(spectrum)),1)));
        tfSingleAmps0(iCh,:,:,iTrig) = specAmp';
        tfSingleITPC0(iCh,:,:,iTrig) = itpc';
    end
end

% blank
for iCh = 1:numel(channels)
    channel = channels(iCh);
    data = squeeze(blankData(:,channel,:))';
    [spectrum,ntaper,freqoi,timeoi] = ft_specest_mtmconvol(data, t/1000, ...
        'timeoi', toi, 'freqoi', foi, 'timwin', t_ftimwin, ...
        'taper', taper, 'dimord', 'chan_time_freqtap');
    specAmp = squeeze(nanmean(abs(spectrum),1)); % mean across trials
    itpc = squeeze(abs(nanmean(exp(1i*angle(spectrum)),1)));
    tfSingleAmps0(iCh,:,:,nTrigs) = specAmp';
    tfSingleITPC0(iCh,:,:,nTrigs) = itpc';
end

% mean across channels
tfSingleAmps1 = squeeze(rd_wmean(tfSingleAmps0,chw,1));
tfSingleITPC = squeeze(rd_wmean(tfSingleITPC0,chw,1));

% normalize by mean amplitude for each frequency
m = nanmean(tfSingleAmps1,2);
tfSingleAmps = tfSingleAmps1./repmat(m,1,numel(timeoi),1)-1; % proportion change from mean

tfSingleAmpsAtt(:,:,1) = nanmean(tfSingleAmps(:,:,plotOrder(1:(nTrigs-1)/2)),3);
tfSingleAmpsAtt(:,:,2) = nanmean(tfSingleAmps(:,:,plotOrder((nTrigs-1)/2+1:end-1)),3);

for iTrig = 1:(nTrigs-1)/2 
    tfSingleAmpsPA(:,:,iTrig) = mean(tfSingleAmps(:,:,iTrig*2-1:iTrig*2),3);
end
t1SinglePADiff = mean(tfSingleAmpsPA(:,:,[1 3]),3)-mean(tfSingleAmpsPA(:,:,[2 4]),3);
t2SinglePADiff = mean(tfSingleAmpsPA(:,:,[1 2]),3)-mean(tfSingleAmpsPA(:,:,[3 4]),3);

% windows around targets
twin = [-600 600];
t1Tidx = find(isneq(timeoi*1000,eventTimes(3)+twin(1))):find(isneq(timeoi*1000,eventTimes(3)+twin(2)));
t2Tidx = find(isneq(timeoi*1000,eventTimes(4)+twin(1))):find(isneq(timeoi*1000,eventTimes(4)+twin(2)));

% calculate pres/abs x att/unattend for each target
stfPAAUT(:,:,1,1) = (tfSingleAmps(:,t1Tidx,1) + tfSingleAmps(:,t1Tidx,5))/2; % present/attended
stfPAAUT(:,:,2,1) = (tfSingleAmps(:,t1Tidx,2) + tfSingleAmps(:,t1Tidx,6))/2; % present/unattended
stfPAAUT(:,:,3,1) = (tfSingleAmps(:,t1Tidx,3) + tfSingleAmps(:,t1Tidx,7))/2; % absent/attended
stfPAAUT(:,:,4,1) = (tfSingleAmps(:,t1Tidx,4) + tfSingleAmps(:,t1Tidx,8))/2; % absent/unattended

stfPAAUT(:,:,1,2) = (tfSingleAmps(:,t2Tidx,2) + tfSingleAmps(:,t2Tidx,4))/2;
stfPAAUT(:,:,2,2) = (tfSingleAmps(:,t2Tidx,1) + tfSingleAmps(:,t2Tidx,3))/2;
stfPAAUT(:,:,3,2) = (tfSingleAmps(:,t2Tidx,6) + tfSingleAmps(:,t2Tidx,8))/2;
stfPAAUT(:,:,4,2) = (tfSingleAmps(:,t2Tidx,5) + tfSingleAmps(:,t2Tidx,7))/2;

% present vs. absent and attended vs. unattended
for iT = 1:2
    stfPAT(:,:,1,iT) = (stfPAAUT(:,:,1,iT) + stfPAAUT(:,:,2,iT))/2; % present
    stfPAT(:,:,2,iT) = (stfPAAUT(:,:,3,iT) + stfPAAUT(:,:,4,iT))/2; % absent

    stfAUT(:,:,1,iT) = (stfPAAUT(:,:,1,iT) + stfPAAUT(:,:,3,iT))/2; % attended
    stfAUT(:,:,2,iT) = (stfPAAUT(:,:,2,iT) + stfPAAUT(:,:,4,iT))/2; % unattended 
end

% combining across T1 and T2
for iPAAU = 1:4
    stfPAAU(:,:,iPAAU) = (stfPAAUT(:,:,iPAAU,1) + stfPAAUT(:,:,iPAAU,2))/2;
end

stfPA(:,:,1) = (stfPAAU(:,:,1) + stfPAAU(:,:,2))/2; % present
stfPA(:,:,2) = (stfPAAU(:,:,3) + stfPAAU(:,:,4))/2; % absent

stfAU(:,:,1) = (stfPAAU(:,:,1) + stfPAAU(:,:,3))/2; % attended
stfAU(:,:,2) = (stfPAAU(:,:,2) + stfPAAU(:,:,4))/2; % unattended

% get values of the time points with respect to target, for plotting
twinvals = toi(t1Tidx)-toi(isneq(toi*1000,eventTimes(3)));


% store results
A.stfTaper = taper;
A.stfFoi = foi;
A.stfTFTimwin = t_ftimwin;
A.stfToi = toi;
A.stfAmps = tfSingleAmps;
A.stfAmpsAtt = tfSingleAmpsAtt;
A.stfAmpsPA = tfSingleAmpsPA;
A.stfPADiff(:,:,1) = t1SinglePADiff;
A.stfPADiff(:,:,2) = t2SinglePADiff;
A.stftwin = twin;
A.stftwinvals = twinvals;
A.stft1Tidx = t1Tidx;
A.stft2Tidx = t2Tidx;
A.stfPAAUT = stfPAAUT;
A.stfPAT = stfPAT;
A.stfAUT = stfAUT;
A.stfPAAU = stfPAAU;
A.stfPA = stfPA;
A.stfAU = stfAU;


% figures
ytick = 10:10:numel(foi);
xtick = 51:xtickint:numel(toi);
clims = [-0.3 0.3]; % [0 70]
diffClims = [-0.2 0.2];
hack = plotOrder;
hack(hack>4) = hack(hack>4)+1;
cmap = flipud(lbmap(64,'RedBlue'));
% cmap = colormap;

if plotFigs
fH = [];
fH(1) = figure;
set(gcf,'Position',tf9FigPos)
for iTrig = 1:nTrigs
    subplot(2,5,hack(iTrig))
    imagesc(tfSingleAmps(:,:,iTrig),clims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    colormap(cmap)
    if iTrig==nTrigs
        xlabel('time (s)')
        ylabel('frequency (Hz)')
    end
    title(trigNames{iTrig})
end
rd_supertitle(['channel' sprintf(' %d', channels) wstrt]);
rd_raiseAxis(gca);

fH(2) = figure;
set(gcf,'Position',tf3FigPos)
attNames = {'attT1','attT2'};
for iAtt = 1:size(tfSingleAmpsAtt,3)
    subplot(1,3,iAtt)
    imagesc(tfSingleAmpsAtt(:,:,iAtt),clims)
    title(attNames{iAtt})
end
subplot(1,3,3)
imagesc(tfSingleAmpsAtt(:,:,2)-tfSingleAmpsAtt(:,:,1),diffClims)
title('attT2 - attT1')
aH = findall(gcf,'type','axes');
for iAx = 1:numel(aH)
    axes(aH(iAx));
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    xlabel('time (s)')
    ylabel('frequency (Hz)')
end
colormap(cmap)
rd_supertitle2(['channel' sprintf(' %d', channels) wstrt]);

fH(3) = figure;
set(gcf,'Position',tf9FigPos)
for iPA = 1:size(tfSingleAmpsPA,3)
    subplot(2,4,iPA)
    imagesc(tfSingleAmpsPA(:,:,iPA),clims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    xlabel('time (s)')
    ylabel('frequency (Hz)')
    title(PANames{iPA})
end
subplot(2,4,5)
imagesc(t1SinglePADiff,diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (s)')
ylabel('frequency (Hz)')
title(sprintf('T1 %s', PADiffNames))
subplot(2,4,6)
imagesc(t2SinglePADiff,diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (s)')
ylabel('frequency (Hz)')
title(sprintf('T2 %s', PADiffNames))
subplot(2,4,7)
imagesc(t2SinglePADiff - t1SinglePADiff,diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (s)')
ylabel('frequency (Hz)')
title(sprintf('T1 vs. T2 %s', PADiffNames))
rd_supertitle2(['channel' sprintf(' %d', channels) wstrt]);
colormap(cmap)

% 9 squares, attended-unattended
fH(4) = figure;
set(gcf,'Position',tf9SquareFigPos)
% T1/T2 x pres/abs
subplot(3,3,1)
imagesc(stfPAAUT(:,:,1,1)-stfPAAUT(:,:,2,1)) % T1-pres-att vs. unatt
ylabel(names{1})
title('T1')
subplot(3,3,2)
imagesc(stfPAAUT(:,:,1,2)-stfPAAUT(:,:,2,2)) % T2-pres-att vs. unatt
title('T2')
subplot(3,3,4)
imagesc(stfPAAUT(:,:,3,1)-stfPAAUT(:,:,4,1)) % T1-abs-att vs. unatt
ylabel(names{2})
subplot(3,3,5)
imagesc(stfPAAUT(:,:,3,2)-stfPAAUT(:,:,4,2)) % T2-abs-att vs. unatt
% ave(T1,T2)
subplot(3,3,3)
imagesc(stfPAAU(:,:,1)-stfPAAU(:,:,2)) % pres-att vs. pres-unatt
title('ave(T1,T2)')
subplot(3,3,6)
imagesc(stfPAAU(:,:,3)-stfPAAU(:,:,4)) % abs-att vs. abs-unatt
% ave(P,A)
subplot(3,3,7)
imagesc(stfAUT(:,:,1,1)-stfAUT(:,:,2,1)) % T1-att vs. T1-unatt 
ylabel(sprintf('ave(%s,\n%s)',names{1},names{2}))
subplot(3,3,8)
imagesc(stfAUT(:,:,1,2)-stfAUT(:,:,2,2)) % T2-att vs. T2-unatt 
% ave(all)
subplot(3,3,9)
imagesc(stfAU(:,:,1)-stfAU(:,:,2)) % att vs. unatt
xlabel('time (s)')
ylabel('frequency (Hz)')
title('ave(all)')
% format subplots
aH = findall(gcf,'type','axes');
paauxtick = [11 61 111];
for iAx = 1:numel(aH)
    axes(aH(iAx));
    rd_timeFreqPlotLabels(twinvals,foi,paauxtick,ytick,0);
    set(gca,'clim',diffClims)
end
colormap(cmap)
rd_supertitle2('attended vs. unattended')

% 9 squares, present-absent
fH(5) = figure;
set(gcf,'Position',tf9SquareFigPos)
% T1/T2 x att/unatt
subplot(3,3,1)
imagesc(stfPAAUT(:,:,1,1)-stfPAAUT(:,:,3,1)) % T1-pres-att vs. abs-att
ylabel('attended')
title('T1')
subplot(3,3,2)
imagesc(stfPAAUT(:,:,1,2)-stfPAAUT(:,:,3,2)) % T2-pres-att vs. abs-att
title('T2')
subplot(3,3,4)
imagesc(stfPAAUT(:,:,2,1)-stfPAAUT(:,:,4,1)) % T1-pres-unatt vs. abs-unatt
ylabel('unattended')
subplot(3,3,5)
imagesc(stfPAAUT(:,:,2,2)-stfPAAUT(:,:,4,2)) % T2-pres-unatt vs. abs-unatt
% ave(T1,T2)
subplot(3,3,3)
imagesc(stfPAAU(:,:,1)-stfPAAU(:,:,3)) % pres-att vs. abs-att
title('ave(T1,T2)')
subplot(3,3,6)
imagesc(stfPAAU(:,:,2)-stfPAAU(:,:,4)) % pres-unatt vs. abs-unatt
% ave(A,U)
subplot(3,3,7)
imagesc(stfPAT(:,:,1,1)-stfPAT(:,:,2,1)) % T1-pres vs. T1-abs 
ylabel('ave(A,U)')
subplot(3,3,8)
imagesc(stfPAT(:,:,1,2)-stfPAT(:,:,2,2)) % T2-pres vs. T2-abs 
% ave(all)
subplot(3,3,9)
imagesc(stfPA(:,:,1)-stfPA(:,:,2)) % pres vs. abs
xlabel('time (s)')
ylabel('frequency (Hz)')
title('ave(all)')
% format subplots
aH = findall(gcf,'type','axes');
paauxtick = [11 61 111];
for iAx = 1:numel(aH)
    axes(aH(iAx));
    rd_timeFreqPlotLabels(twinvals,foi,paauxtick,ytick,0);
    set(gca,'clim',diffClims)
end
colormap(cmap)
rd_supertitle2(sprintf('%s vs. %s',names{1},names{2}))

% 6 squares, present
fH(6) = figure;
set(gcf,'Position',tf6SquareFigPos)
% T1/T2 x att/unatt
subplot(2,3,1)
imagesc(stfPAAUT(:,:,1,1)) % T1-pres-att
ylabel('attended')
title('T1')
subplot(2,3,2)
imagesc(stfPAAUT(:,:,1,2)) % T2-pres-att
title('T2')
subplot(2,3,4)
imagesc(stfPAAUT(:,:,2,1)) % T1-pres-unatt
ylabel('unattended')
subplot(2,3,5)
imagesc(stfPAAUT(:,:,2,2)) % T2-pres-unatt
% ave(T1,T2)
subplot(2,3,3)
imagesc(stfPAAU(:,:,1)) % pres-att
title('ave(T1,T2)')
subplot(2,3,6)
imagesc(stfPAAU(:,:,2)) % pres-unatt
xlabel('time (s)')
ylabel('frequency (Hz)')
title('ave(all)')
% format subplots
aH = findall(gcf,'type','axes');
paauxtick = [11 61 111];
for iAx = 1:numel(aH)
    axes(aH(iAx));
    rd_timeFreqPlotLabels(twinvals,foi,paauxtick,ytick,0);
    set(gca,'clim',diffClims)
end
colormap(cmap)
rd_supertitle2(names{1})

% 6 squares, absent
fH(7) = figure;
set(gcf,'Position',tf6SquareFigPos)
% T1/T2 x att/unatt
subplot(2,3,1)
imagesc(stfPAAUT(:,:,3,1)) % T1-abs-att
ylabel('attended')
title('T1')
subplot(2,3,2)
imagesc(stfPAAUT(:,:,3,2)) % T2-abs-att
title('T2')
subplot(2,3,4)
imagesc(stfPAAUT(:,:,4,1)) % T1-abs-unatt
ylabel('unattended')
subplot(2,3,5)
imagesc(stfPAAUT(:,:,4,2)) % T2-abs-unatt
% ave(T1,T2)
subplot(2,3,3)
imagesc(stfPAAU(:,:,3)) % abs-att
title('ave(T1,T2)')
subplot(2,3,6)
imagesc(stfPAAU(:,:,4)) % abs-unatt
xlabel('time (s)')
ylabel('frequency (Hz)')
title('ave(all)')
% format subplots
aH = findall(gcf,'type','axes');
paauxtick = [11 61 111];
for iAx = 1:numel(aH)
    axes(aH(iAx));
    rd_timeFreqPlotLabels(twinvals,foi,paauxtick,ytick,0);
    set(gca,'clim',diffClims)
end
colormap(cmap)
rd_supertitle2(names{2})
end

if saveFigs
    if numel(channels)==1
        figPrefix = sprintf('im_ch%d', channels);
    else
        figPrefix = ['im_ch' sprintf('%d_', channels(1:end-1)) sprintf('%d', channels(end)) wstr];
    end
    rd_saveAllFigs(fH, {'timeFreqSingleByCond','timeFreqSingleAtt','timeFreqSinglePA','timeFreqSingleAUDiff','timeFreqSinglePADiff','timeFreqSinglePresent','timeFreqSingleAbsent'}, figPrefix, figDir)
end

%% inter-trial phase coherence
tfSingleITPCAtt(:,:,1) = nanmean(tfSingleITPC(:,:,plotOrder(1:(nTrigs-1)/2)),3);
tfSingleITPCAtt(:,:,2) = nanmean(tfSingleITPC(:,:,plotOrder((nTrigs-1)/2+1:end-1)),3);

for iTrig = 1:(nTrigs-1)/2 
    tfSingleITPCPA(:,:,iTrig) = mean(tfSingleITPC(:,:,iTrig*2-1:iTrig*2),3);
end
t1SingleITPCPADiff = mean(tfSingleITPCPA(:,:,[1 3]),3)-mean(tfSingleITPCPA(:,:,[2 4]),3);
t2SingleITPCPADiff = mean(tfSingleITPCPA(:,:,[1 2]),3)-mean(tfSingleITPCPA(:,:,[3 4]),3);

% in twin
% calculate pres/abs x att/unattend for each target
itpcPAAUT(:,:,1,1) = (tfSingleITPC(:,t1Tidx,1) + tfSingleITPC(:,t1Tidx,5))/2; % present/attended
itpcPAAUT(:,:,2,1) = (tfSingleITPC(:,t1Tidx,2) + tfSingleITPC(:,t1Tidx,6))/2; % present/unattended
itpcPAAUT(:,:,3,1) = (tfSingleITPC(:,t1Tidx,3) + tfSingleITPC(:,t1Tidx,7))/2; % absent/attended
itpcPAAUT(:,:,4,1) = (tfSingleITPC(:,t1Tidx,4) + tfSingleITPC(:,t1Tidx,8))/2; % absent/unattended

itpcPAAUT(:,:,1,2) = (tfSingleITPC(:,t2Tidx,2) + tfSingleITPC(:,t2Tidx,4))/2;
itpcPAAUT(:,:,2,2) = (tfSingleITPC(:,t2Tidx,1) + tfSingleITPC(:,t2Tidx,3))/2;
itpcPAAUT(:,:,3,2) = (tfSingleITPC(:,t2Tidx,6) + tfSingleITPC(:,t2Tidx,8))/2;
itpcPAAUT(:,:,4,2) = (tfSingleITPC(:,t2Tidx,5) + tfSingleITPC(:,t2Tidx,7))/2;

% present vs. absent and attended vs. unattended
for iT = 1:2
    itpcPAT(:,:,1,iT) = (itpcPAAUT(:,:,1,iT) + itpcPAAUT(:,:,2,iT))/2; % present
    itpcPAT(:,:,2,iT) = (itpcPAAUT(:,:,3,iT) + itpcPAAUT(:,:,4,iT))/2; % absent

    itpcAUT(:,:,1,iT) = (itpcPAAUT(:,:,1,iT) + itpcPAAUT(:,:,3,iT))/2; % attended
    itpcAUT(:,:,2,iT) = (itpcPAAUT(:,:,2,iT) + itpcPAAUT(:,:,4,iT))/2; % unattended 
end

% combining across T1 and T2
for iPAAU = 1:4
    itpcPAAU(:,:,iPAAU) = (itpcPAAUT(:,:,iPAAU,1) + itpcPAAUT(:,:,iPAAU,2))/2;
end

itpcPA(:,:,1) = (itpcPAAU(:,:,1) + itpcPAAU(:,:,2))/2; % present
itpcPA(:,:,2) = (itpcPAAU(:,:,3) + itpcPAAU(:,:,4))/2; % absent

itpcAU(:,:,1) = (itpcPAAU(:,:,1) + itpcPAAU(:,:,3))/2; % attended
itpcAU(:,:,2) = (itpcPAAU(:,:,2) + itpcPAAU(:,:,4))/2; % unattended

A.stfITPCAmps = tfSingleITPC;
A.stfITPCAtt = tfSingleITPCAtt;
A.stfITPCPA = tfSingleITPCPA;
A.stfITPCPADiff(:,:,1) = t1SingleITPCPADiff;
A.stfITPCPADiff(:,:,2) = t2SingleITPCPADiff;
A.stfITPCPAAUT = itpcPAAUT;
A.stfITPCPAT = itpcPAT;
A.stfITPCAUT = itpcAUT;
A.stfITPCPAAU = itpcPAAU;
A.stfITPCPA = itpcPA;
A.stfITPCAU = itpcAU;

% figures
clims = [0 1]; % [0 70]
diffClims = [-0.2 0.2];
cmap = parula;

if plotFigs
fH = [];
fH(1) = figure;
set(gcf,'Position',tf9FigPos)
for iTrig = 1:nTrigs
    subplot(2,5,hack(iTrig))
    imagesc(tfSingleITPC(:,:,iTrig),clims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    colormap(cmap)
    if iTrig==nTrigs
        xlabel('time (s)')
        ylabel('frequency (Hz)')
    end
    title(trigNames{iTrig})
end
rd_supertitle(['channel' sprintf(' %d', channels) wstrt]);
rd_raiseAxis(gca);

fH(2) = figure;
set(gcf,'Position',tf3FigPos)
attNames = {'attT1','attT2'};
for iAtt = 1:size(tfSingleITPCAtt,3)
    subplot(1,3,iAtt)
    imagesc(tfSingleITPCAtt(:,:,iAtt),clims)
    title(attNames{iAtt})
end
subplot(1,3,3)
imagesc(tfSingleITPCAtt(:,:,2)-tfSingleITPCAtt(:,:,1))%,diffClims)
title('attT2 - attT1')
aH = findall(gcf,'type','axes');
for iAx = 1:numel(aH)
    axes(aH(iAx));
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    xlabel('time (s)')
    ylabel('frequency (Hz)')
end
colormap(cmap)
rd_supertitle2(['channel' sprintf(' %d', channels) wstrt]);
end

if saveFigs
    if numel(channels)==1
        figPrefix = sprintf('im_ch%d', channels);
    else
        figPrefix = ['im_ch' sprintf('%d_', channels(1:end-1)) sprintf('%d', channels(end)) wstr];
    end
    rd_saveAllFigs(fH, {'timeFreqITPCByCond','timeFreqITPCAtt'}, figPrefix, figDir)
end

%% save analysis
if saveAnalysis
    save(analysisFileName, 'A')
end

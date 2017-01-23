function rd_TADetectDiscrimDVA(exptDir, sessionDir, fileBase, analStr, ~, trialSelection, respTargetSelection)

% stability - dva
% whole brain single trial analysis (first part as in SSVEF5)

%% Setup
if nargin==0 || ~exist('exptDir','var')
    exptDir = '/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG';
    sessionDir = 'R0817_20150504';
    fileBase = 'R0817_TADeDi_5.4.15';
    analStr = 'ebi'; % '', 'ebi', etc.
    trialSelection = 'all'; % 'all','validCorrect'
    respTargetSelection = ''; % '','T1Resp','T2Resp'
end

channelSelectionStr = 'wholebrain';

dataDir = sprintf('%s/%s', exptDir, sessionDir);
matDir = sprintf('%s/mat', dataDir);

switch analStr
    case ''
        savename = sprintf('%s/%s_ssvef_workspace.mat', matDir, fileBase);
        analysisFileName = sprintf('%s/analysis_singleTrials_%s_%s_%sTrials%s_dva.mat', matDir, fileBase, channelSelectionStr, trialSelection, respTargetSelection);
    otherwise
        savename = sprintf('%s/%s_%s_ssvef_workspace.mat', matDir, fileBase, analStr);
        analysisFileName = sprintf('%s/analysis_singleTrials_%s_%s_%s_%sTrials%s_dva.mat', matDir, fileBase, analStr, channelSelectionStr, trialSelection, respTargetSelection);
end

%% Get the data
load(savename)

%% Update behav
behav = behavior(behav);

%% Settings after loading the data
saveAnalysis = 1;
saveFigs = 0;
plotFigs = 1;

excludeTrialsFt = 1;
excludeSaturatedEpochs = 0;

channels = 1:157;

%% Plotting setup
set(0,'defaultLineLineWidth',1)

%% Store settings for this analysis
A.fileBase = fileBase;
A.analStr = analStr;
A.trialSelection = trialSelection;
A.respTargetSelection = respTargetSelection;
A.excludeTrialsFt = excludeTrialsFt;
A.excludeSaturatedEpochs = excludeSaturatedEpochs;
A.channels = channels;
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
            analysisFileName = sprintf('%s/analysis_singleTrials_%s_ft_%s_%sTrials%s_dva.mat', matDir, fileBase, channelSelectionStr, trialSelection, respTargetSelection);
        otherwise
            analysisFileName = sprintf('%s/analysis_singleTrials_%s_%s_ft_%s_%sTrials%s_dva.mat', matDir, fileBase, analStr, channelSelectionStr, trialSelection, respTargetSelection);
    end
end

%% Exclude channels manually rejected with ft
% NEED TO ADD

%% Make figDir if needed
figDir = sprintf('%s_singleTrials_%s_%sTrials%s', figDir, channelSelectionStr, trialSelection, respTargetSelection);

if ~exist(figDir,'dir') && saveFigs
    mkdir(figDir)
end

%% Trial selection
cueCondIdx = strcmp(behav.responseData_labels, 'cue condition');
nTrials = size(behav.responseData_all,1);

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

blankCond = 1;
wBlank = behav.responseData_all(:,cueCondIdx) == blankCond;

%% remove unselected and blank trials
trigDataSelected = trigData;
trigDataSelected(:,:,~wSelect) = NaN;
trigDataSelected(:,:,wBlank) = NaN;

excludedTrials = isnan(squeeze(trigDataSelected(1,1,:)));
data = trigDataSelected(:,:,~excludedTrials);

%% visualize one trial
if plotFigs
    figure('Position',[200 500 1000 400])
    imagesc(data(:,:,1)')
    for iEv = 1:numel(eventTimes)
        vline(find(t==eventTimes(iEv)),'k');
    end
    xlabel('time (ms)')
    ylabel('channel')
end

%% calculate R, dva, and L2 norm for each trial
% initialize
R = nan(size(data,1),size(data,3));
L2 = nan(size(R));
n = 100; % n time points in each bin

for iTrial = 1:size(data,3)
    if mod(iTrial,10)==0
    fprintf('trial %d\n', iTrial)
    end
    for idx = n/2+1:size(data,1)-n/2; % center of each bin
        d = data(idx-n/2:idx+n/2-1,:,iTrial)'; % n time points, all channels, one trial
        % normalize each vector
        for iTime = 1:size(d,2)
            dnorm(:,iTime) = d(:,iTime)/norm(d(:,iTime));
        end
        R(idx,iTrial) = norm(sum(dnorm,2))/n;
        L2(idx,iTrial) = norm(data(idx,:,iTrial));
    end
end
dva = 1-R;

A.R = R;
A.dva = dva;
A.L2 = L2;

%% plot
if plotFigs
    figure
    hold on
    plot(mean(dva,2))
    for iEv = 1:numel(eventTimes)
        vline(find(t==eventTimes(iEv)),'k');
    end
    xlabel('time (ms)')
    ylabel('dva')
    
    figure
    hold on
    plot(mean(L2,2))
    for iEv = 1:numel(eventTimes)
        vline(find(t==eventTimes(iEv)),'k');
    end
    xlabel('time (ms)')
    ylabel('L2 norm')
end

%% save analysis
if saveAnalysis
    save(analysisFileName, 'A') % '-v7.3'
end



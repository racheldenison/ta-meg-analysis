function rd_TADetectDiscrimSSVEF5(exptDir, sessionDir, fileBase, analStr, ssvefFreq, nTopChannels, iqrThresh, weightChannels, trialSelection, respTargetSelection)

% whole brain single trial analysis

%% Setup
if nargin==0 || ~exist('exptDir','var')
    exptDir = '/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG';
    sessionDir = 'R0817_20150504';
    fileBase = 'R0817_TADeDi_5.4.15';
    analStr = 'ebi'; % '', 'ebi', etc.
    ssvefFreq = 30;
    trialSelection = 'all'; % 'all','validCorrect'
    respTargetSelection = ''; % '','T1Resp','T2Resp'
end

channelSelectionStr = 'wholebrain';

dataDir = sprintf('%s/%s', exptDir, sessionDir);
matDir = sprintf('%s/mat', dataDir);

switch analStr
    case ''
        savename = sprintf('%s/%s_ssvef_workspace.mat', matDir, fileBase);
        analysisFileName = sprintf('%s/analysis_singleTrials_%s_%s_%sTrials%s_%dHz.mat', matDir, fileBase, channelSelectionStr, trialSelection, respTargetSelection, ssvefFreq);
    otherwise
        savename = sprintf('%s/%s_%s_ssvef_workspace.mat', matDir, fileBase, analStr);
        analysisFileName = sprintf('%s/analysis_singleTrials_%s_%s_%s_%sTrials%s_%dHz.mat', matDir, fileBase, analStr, channelSelectionStr, trialSelection, respTargetSelection, ssvefFreq);
end

% load data header for plotting topologies
load data/data_hdr.mat

%% Get the data
load(savename)

%% Update behav
behav = behavior(behav);

%% Settings after loading the data
saveAnalysis = 0;
saveFigs = 0;
plotFigs = 1;

excludeTrialsFt = 1;
excludeSaturatedEpochs = 0;

channels = 1:157;

%% Store settings for this analysis
A.fileBase = fileBase;
A.analStr = analStr;
A.excludeTrialsFt = excludeTrialsFt;
A.excludeSaturatedEpochs = excludeSaturatedEpochs;
A.ssvefFreq = ssvefFreq;
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
            analysisFileName = sprintf('%s/analysis_singleTrials_%s_ft_%s_%sTrials%s_%dHz.mat', matDir, fileBase, channelSelectionStr, trialSelection, respTargetSelection, ssvefFreq);
        otherwise
            analysisFileName = sprintf('%s/analysis_singleTrials_%s_%s_ft_%s_%sTrials%s_%dHz.mat', matDir, fileBase, analStr, channelSelectionStr, trialSelection, respTargetSelection, ssvefFreq);
    end
end

%% Exclude channels manually rejected with ft
% NEED TO ADD

%% Make figDir if needed
figDir = sprintf('%s_singleTrials_%s_%sTrials%s', figDir, channelSelectionStr, trialSelection, respTargetSelection);

if ~exist(figDir,'dir') && saveFigs
    mkdir(figDir)
end

%% Organize trials into conditions
cueCondIdx = strcmp(behav.responseData_labels, 'cue condition');
t1CondIdx = strcmp(behav.responseData_labels, 'target type T1');
t2CondIdx = strcmp(behav.responseData_labels, 'target type T2');
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

trigDataSelected = trigData; % make a copy so we use it for condData but not blankData
trigDataSelected(:,:,wSelect~=1)=NaN;

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

% mean across trials
trigMeanMean = squeeze(nanmean(trigMean,3));

A.nTrialsCond = nTrialsCond;
A.trigMean = trigMean;
A.trigMeanMean = trigMeanMean;

%% Wavelet
switch ssvefFreq
    case 30
        width = 12; % 12 for 30 Hz, 16 for 40 Hz gives 127 ms duration, 5 Hz bandwidth
    case 40
        width = 16;
    otherwise
        error('ssvefFreq not recognized')
end

waveletOption = 'trialmean'; % 'trialmean', 'singletrials'
switch waveletOption
    case 'trialmean'
        wAmps = [];
        foi = ssvefFreq;
        for iTrig = 1:nTrigs
            data = trigMeanMean(:,:,iTrig)'; % channels by samples
            [spectrum,freqoi,timeoi] = ft_specest_wavelet(data, t/1000, 'freqoi', foi, 'width', width);
            specAmp = abs(squeeze(spectrum));
            
            if all(size(specAmp)>1) % if two-dimensional
                wAmp = specAmp;
            else
                wAmp = specAmp';
            end
            wAmps(:,:,iTrig) = wAmp';
        end
    case 'singletrials'
        wAmps0 = [];
        foi = ssvefFreq;
        for iTrig = 1:nTrigs
            for iTrial = 1:nTrialsPerCond
                data = trigMean(:,:,iTrial,iTrig)'; % channels by samples
                [spectrum,freqoi,timeoi] = ft_specest_wavelet(data, t/1000, 'freqoi', foi, 'width', width);
                specAmp = abs(squeeze(spectrum));
                
                if all(size(specAmp)>1) % if two-dimensional
                    wAmp = specAmp;
                else
                    wAmp = specAmp';
                end
                wAmps0(:,:,iTrial,iTrig) = wAmp';
            end
        end
        wAmps = squeeze(nanmean(wAmps0,3)); % mean across trials
    otherwise
        error('waveletOption not recognized')
end

%% PAAU
twin = [-600 600];
t1Tidx = find(t==eventTimes(3)+twin(1)):find(t==eventTimes(3)+twin(2));
t2Tidx = find(t==eventTimes(4)+twin(1)):find(t==eventTimes(4)+twin(2));

% calculate pres/abs x att/unattend for each target, groupData
wPAAUT(:,:,1,1) = mean(wAmps(t1Tidx,:,[1 5]),3); % present/attended
wPAAUT(:,:,2,1) = mean(wAmps(t1Tidx,:,[2 6]),3); % present/unattended
wPAAUT(:,:,3,1) = mean(wAmps(t1Tidx,:,[3 7]),3); % absent/attended
wPAAUT(:,:,4,1) = mean(wAmps(t1Tidx,:,[4 8]),3); % absent/unattended

wPAAUT(:,:,1,2) = mean(wAmps(t2Tidx,:,[2 4]),3);
wPAAUT(:,:,2,2) = mean(wAmps(t2Tidx,:,[1 3]),3);
wPAAUT(:,:,3,2) = mean(wAmps(t2Tidx,:,[6 8]),3);
wPAAUT(:,:,4,2) = mean(wAmps(t2Tidx,:,[5 7]),3);

% present vs. absent and attended vs. unattended
for iT = 1:2
    wPAT(:,:,1,iT) = mean(wPAAUT(:,:,[1 2],iT),3); % present
    wPAT(:,:,2,iT) = mean(wPAAUT(:,:,[3 4],iT),3); % absent

    wAUT(:,:,1,iT) = mean(wPAAUT(:,:,[1 3],iT),3); % attended
    wAUT(:,:,2,iT) = mean(wPAAUT(:,:,[2 4],iT),3); % unattended 
end

% combining across T1 and T2
for iPAAU = 1:4
    wPAAU(:,:,iPAAU) = mean(wPAAUT(:,:,iPAAU,[1 2]),4);
end

wPA(:,:,1) = mean(wPAAU(:,:,[1 2]),3); % present
wPA(:,:,2) = mean(wPAAU(:,:,[3 4]),3); % absent

wAU(:,:,1) = mean(wPAAU(:,:,[1 3]),3); % attended
wAU(:,:,2) = mean(wPAAU(:,:,[2 4]),3); % unattended

% finally, ampsAll
wAmpsAll = squeeze(mean(wAmps(:,:,1:end-1),3));

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
A.wAmpsAll = wAmpsAll;

%% Plot wavelet topo
paauNames = {'P-att','P-unatt','A-att','A-unatt'};
auNames = {'att','unatt'};
paNames = {'P','A'};

%% movie
twindow = twin(1):twin(end);
tstep = 10;
iPAAU = 2;
iT = 1;
figure
for iTime = 1:tstep:numel(twindow)
    fH = ssm_plotOnMesh(wPAAUT(iTime,:,iPAAU,iT), ...
        sprintf('wPAAUT t=%d',twindow(iTime)),[], data_hdr, '2d');
    set(gca,'CLim',[0 200])
    pause(.1)
end

%% time bins
% PAAU
% figPos = [32 250 2400 450];
figPos = [32 150 2400 750];
clims = [0 200];
nBins = 12;
binSize = round(numel(twindow)/nBins);
for iT = 1:2
    figure('Position',figPos);
    for iPAAU = 1:4
        for iBin = 1:nBins
            subplot(4,nBins,iBin + nBins*(iPAAU-1))
            tidx = (1:binSize+1) + (iBin-1)*binSize;
            str = sprintf('wPAAUT, %s, t=[%d %d]',paauNames{iPAAU}, twindow(tidx(1)), twindow(tidx(end)));
            vals = mean(wPAAUT(tidx,:,iPAAU,iT),1);
            fH = ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
            set(gca,'CLim',clims)
        end
    end
    rd_supertitle2(sprintf('T%d',iT))
end

% AU
figPos = [32 250 2400 450];
% figPos = [35 150 2400 750];
clims = [0 200];
nBins = 12;
binSize = round(numel(twindow)/nBins);
for iT = 1:2
    figure('Position',figPos);
    for iAU = 1:2
        for iBin = 1:nBins
            subplot(2,nBins,iBin + nBins*(iAU-1))
            tidx = (1:binSize+1) + (iBin-1)*binSize;
            str = sprintf('wAUT, %s, t=[%d %d]',auNames{iAU}, twindow(tidx(1)), twindow(tidx(end)));
            vals = mean(wAUT(tidx,:,iAU,iT),1);
            fH = ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
            set(gca,'CLim',clims)
        end
    end
    rd_supertitle2(sprintf('T%d',iT))
end

% PA
figPos = [32 250 2400 450];
% figPos = [35 150 2400 750];
clims = [0 200];
nBins = 12;
binSize = round(numel(twindow)/nBins);
for iT = 1:2
    figure('Position',figPos);
    for iPA = 1:2
        for iBin = 1:nBins
            subplot(2,nBins,iBin + nBins*(iPA-1))
            tidx = (1:binSize+1) + (iBin-1)*binSize;
            str = sprintf('wPAT, %s, t=[%d %d]',paNames{iPA}, twindow(tidx(1)), twindow(tidx(end)));
            vals = mean(wPAT(tidx,:,iPA,iT),1);
            fH = ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
            set(gca,'CLim',clims)
        end
    end
    rd_supertitle2(sprintf('T%d',iT))
end

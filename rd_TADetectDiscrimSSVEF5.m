function rd_TADetectDiscrimSSVEF5(exptDir, sessionDir, fileBase, analStr, ssvefFreq, trialSelection, respTargetSelection)

% whole brain single trial analysis

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
            trialSelection = 'all'; % 'all','validCorrect', etc
            respTargetSelection = ''; % '','T1Resp','T2Resp'
            
        case 'TAContrast'
            exptDir = '/Local/Users/denison/Data/TAContrast/MEG';
            sessionDir = 'R0817_20171019';
            fileBase = 'R0817_TACont_10.19.17';
            analStr = 'ebi'; % '', 'ebi', etc.
            ssvefFreq = 20;
            trialSelection = 'all'; % 'all','validCorrect', etc
            respTargetSelection = ''; % '','T1Resp','T2Resp'
            
        case 'TANoise'
            exptDir = '/Local/Users/denison/Data/TANoise/MEG';
            sessionDir = 'R0817_20171212';
            fileBase = 'R0817_TANoise_12.12.17';
            analStr = 'ebi'; % '', 'ebi', etc.
            ssvefFreq = 20;
            trialSelection = 'all'; % 'all','validCorrect', etc
            respTargetSelection = ''; % '','T1Resp','T2Resp'
            
        otherwise
            error('exptType not recognized')
    end
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
cfg = [];
layout = ft_prepare_layout(cfg, data_hdr);

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

channels = 1:157;

%% Plotting setup
plotOrder = [1 5 3 7 2 6 4 8 9];

tf9FigPos = [0 250 1280 580];
tf3FigPos = [200 475 980 330];
paau3FigPos = [800 30 450 830];
paau6FigPos = [0 90 1000 800];
tf9SquareFigPos = [50 50 850 850];
tf6SquareFigPos = [50 50 850 530];

set(0,'defaultLineLineWidth',1)

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

% mean across trials
trigMeanMean = squeeze(nanmean(trigMean,3));

A.nTrialsCond = nTrialsCond;
% A.trigMean = trigMean; % very big, so not saving
A.trigMeanMean = trigMeanMean;

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

A.waveletOption = waveletOption;
A.wAmps = wAmps;

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
switch exptType
    case 'TADetectDiscrim'
        names = {'target present','target absent'};
        paauNames = {'P-att','P-unatt','A-att','A-unatt'};
        paNames = {'T1p-T2p','T1a-T2p','T1p-T2a','T1a-T2a'};
        paDiffNames = 'P-A';
        xtickint = 50;
    case 'TAContrast'
        names = {'target decrement','target increment'};
        paauNames = {'D-att','D-unatt','I-att','I-unatt'};
        paNames = {'T1d-T2d','T1i-T2d','T1d-T2i','T1i-T2i'};
        paDiffNames = 'D-I';
        xtickint = 100;
end

auNames = {'att','unatt'};
clims = [0 800];
diffClims = [-200 200];
twindow = twin(1):twin(end);
nBins = 6;
binSize = round(numel(twindow)/nBins);
load parula
cmap = flipud(lbmap(64,'RedBlue'));
% set(groot,'DefaultFigureColormap',cmap)

%% movie
% if plotFigs
% tstep = 10;
% iPAAU = 1;
% iT = 1;
% figure
% for iTime = 1:tstep:numel(twindow)
%     fH = ssm_plotOnMesh(wPAAUT(iTime,:,iPAAU,iT), ...
%         sprintf('wPAAUT t=%d',twindow(iTime)),[], data_hdr, '2d');
%     set(gca,'CLim',clims)
%     pause(.1)
% end
% end

%% time bins
if plotFigs
fH = [];
% PAAU
figPos = [32 150 nBins*200 750];
for iT = 1:2
    fH(0+iT) = figure('Position',figPos);
    for iPAAU = 1:4
        for iBin = 1:nBins
            subplot(4,nBins,iBin + nBins*(iPAAU-1))
            tidx = (1:binSize+1) + (iBin-1)*binSize;
            str = sprintf('wPAAUT, %s, t=[%d %d]',paauNames{iPAAU}, twindow(tidx(1)), twindow(tidx(end)));
            vals = mean(wPAAUT(tidx,:,iPAAU,iT),1);
            ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
            set(gca,'CLim',clims)
        end
    end
    colormap(parula)
    rd_supertitle2(sprintf('T%d',iT))
end

% AU
figPos = [32 250 nBins*200 650];
for iT = 1:2
    fH(2+iT) = figure('Position',figPos);
    for iAU = 1:2
        for iBin = 1:nBins
            subplot(3,nBins,iBin + nBins*(iAU-1))
            tidx = (1:binSize+1) + (iBin-1)*binSize;
            str = sprintf('wAUT, %s, t=[%d %d]',auNames{iAU}, twindow(tidx(1)), twindow(tidx(end)));
            vals = mean(wAUT(tidx,:,iAU,iT),1);
            ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
            set(gca,'CLim',clims)
            colormap(parula)
            freezeColors
        end
    end
    for iBin = 1:nBins
        subplot(3,nBins,iBin + nBins*2)
        tidx = (1:binSize+1) + (iBin-1)*binSize;
        str = sprintf('wAUT, A-U, t=[%d %d]',twindow(tidx(1)), twindow(tidx(end)));
        vals = mean((wAUT(tidx,:,1,iT) - wAUT(tidx,:,2,iT)),1);
        ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
        set(gca,'CLim',diffClims)
    end
    colormap(cmap)
    rd_supertitle2(sprintf('T%d',iT))
end

% PA
figPos = [32 250 nBins*200 650];
for iT = 1:2
    fH(4+iT) = figure('Position',figPos);
    for iPA = 1:2
        for iBin = 1:nBins
            subplot(3,nBins,iBin + nBins*(iPA-1))
            tidx = (1:binSize+1) + (iBin-1)*binSize;
            str = sprintf('wPAT, %s, t=[%d %d]',paNames{iPA}, twindow(tidx(1)), twindow(tidx(end)));
            vals = mean(wPAT(tidx,:,iPA,iT),1);
            ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
            set(gca,'CLim',clims)
            colormap(parula)
            freezeColors
        end
    end
    for iBin = 1:nBins
        subplot(3,nBins,iBin + nBins*2)
        tidx = (1:binSize+1) + (iBin-1)*binSize;
        str = sprintf('wPAT, %s, t=[%d %d]',paDiffNames, twindow(tidx(1)), twindow(tidx(end)));
        vals = mean((wPAT(tidx,:,1,iT) - wPAT(tidx,:,2,iT)),1);
        ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
        set(gca,'CLim',diffClims)
    end
    colormap(cmap)
    rd_supertitle2(sprintf('T%d',iT))
end
end

if saveFigs
    figPrefix = ['map_wholebrain' sprintf('%dHz', ssvefFreq)];
    rd_saveAllFigs(fH, {'wPAAUT1','wPAAUT2','wAUT1','wAUT2','wPAT1','wPAT2'}, figPrefix, figDir)
end

%% Time-frequency - single trials
taper          = 'hanning';
foi            = 1:50;
t_ftimwin      = 10 ./ foi;
toi            = tstart/1000:0.01:tstop/1000;
tfSingleAmps0 = [];
tfSingleITPC = [];
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
        tfSingleITPC(iCh,:,:,iTrig) = itpc';
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
    tfSingleITPC(iCh,:,:,nTrigs) = itpc';
end

% normalize by mean amplitude for each frequency **in blank condition**
m = nanmean(tfSingleAmps0(:,:,:,end),3);
tfSingleAmps = tfSingleAmps0./repmat(m,1,1,numel(timeoi),nTrigs)-1; % proportion change from mean

tfSingleAmpsAtt(:,:,:,1) = nanmean(tfSingleAmps(:,:,:,plotOrder(1:(nTrigs-1)/2)),4);
tfSingleAmpsAtt(:,:,:,2) = nanmean(tfSingleAmps(:,:,:,plotOrder((nTrigs-1)/2+1:end-1)),4);

for iTrig = 1:(nTrigs-1)/2 
    tfSingleAmpsPA(:,:,:,iTrig) = mean(tfSingleAmps(:,:,:,iTrig*2-1:iTrig*2),4);
end
t1SinglePADiff = mean(tfSingleAmpsPA(:,:,:,[1 3]),4)-mean(tfSingleAmpsPA(:,:,:,[2 4]),4);
t2SinglePADiff = mean(tfSingleAmpsPA(:,:,:,[1 2]),4)-mean(tfSingleAmpsPA(:,:,:,[3 4]),4);

% windows around targets
twin = [-600 600];
t1Tidx = find(isneq(timeoi*1000,eventTimes(3)+twin(1))):find(isneq(timeoi*1000,eventTimes(3)+twin(2)));
t2Tidx = find(isneq(timeoi*1000,eventTimes(4)+twin(1))):find(isneq(timeoi*1000,eventTimes(4)+twin(2)));

% calculate pres/abs x att/unattend for each target
stfPAAUT(:,:,:,1,1) = (tfSingleAmps(:,:,t1Tidx,1) + tfSingleAmps(:,:,t1Tidx,5))/2; % present/attended
stfPAAUT(:,:,:,2,1) = (tfSingleAmps(:,:,t1Tidx,2) + tfSingleAmps(:,:,t1Tidx,6))/2; % present/unattended
stfPAAUT(:,:,:,3,1) = (tfSingleAmps(:,:,t1Tidx,3) + tfSingleAmps(:,:,t1Tidx,7))/2; % absent/attended
stfPAAUT(:,:,:,4,1) = (tfSingleAmps(:,:,t1Tidx,4) + tfSingleAmps(:,:,t1Tidx,8))/2; % absent/unattended

stfPAAUT(:,:,:,1,2) = (tfSingleAmps(:,:,t2Tidx,2) + tfSingleAmps(:,:,t2Tidx,4))/2;
stfPAAUT(:,:,:,2,2) = (tfSingleAmps(:,:,t2Tidx,1) + tfSingleAmps(:,:,t2Tidx,3))/2;
stfPAAUT(:,:,:,3,2) = (tfSingleAmps(:,:,t2Tidx,6) + tfSingleAmps(:,:,t2Tidx,8))/2;
stfPAAUT(:,:,:,4,2) = (tfSingleAmps(:,:,t2Tidx,5) + tfSingleAmps(:,:,t2Tidx,7))/2;

% present vs. absent and attended vs. unattended
for iT = 1:2
    stfPAT(:,:,:,1,iT) = (stfPAAUT(:,:,:,1,iT) + stfPAAUT(:,:,:,2,iT))/2; % present
    stfPAT(:,:,:,2,iT) = (stfPAAUT(:,:,:,3,iT) + stfPAAUT(:,:,:,4,iT))/2; % absent

    stfAUT(:,:,:,1,iT) = (stfPAAUT(:,:,:,1,iT) + stfPAAUT(:,:,:,3,iT))/2; % attended
    stfAUT(:,:,:,2,iT) = (stfPAAUT(:,:,:,2,iT) + stfPAAUT(:,:,:,4,iT))/2; % unattended 
end

% combining across T1 and T2
for iPAAU = 1:4
    stfPAAU(:,:,:,iPAAU) = (stfPAAUT(:,:,:,iPAAU,1) + stfPAAUT(:,:,:,iPAAU,2))/2;
end

stfPA(:,:,:,1) = (stfPAAU(:,:,:,1) + stfPAAU(:,:,:,2))/2; % present
stfPA(:,:,:,2) = (stfPAAU(:,:,:,3) + stfPAAU(:,:,:,4))/2; % absent

stfAU(:,:,:,1) = (stfPAAU(:,:,:,1) + stfPAAU(:,:,:,3))/2; % attended
stfAU(:,:,:,2) = (stfPAAU(:,:,:,2) + stfPAAU(:,:,:,4))/2; % unattended

% inter-trial phase coherence
% calculate pres/abs x att/unattend for each target
itpcPAAUT(:,:,:,1,1) = (tfSingleITPC(:,:,t1Tidx,1) + tfSingleITPC(:,:,t1Tidx,5))/2; % present/attended
itpcPAAUT(:,:,:,2,1) = (tfSingleITPC(:,:,t1Tidx,2) + tfSingleITPC(:,:,t1Tidx,6))/2; % present/unattended
itpcPAAUT(:,:,:,3,1) = (tfSingleITPC(:,:,t1Tidx,3) + tfSingleITPC(:,:,t1Tidx,7))/2; % absent/attended
itpcPAAUT(:,:,:,4,1) = (tfSingleITPC(:,:,t1Tidx,4) + tfSingleITPC(:,:,t1Tidx,8))/2; % absent/unattended

itpcPAAUT(:,:,:,1,2) = (tfSingleITPC(:,:,t2Tidx,2) + tfSingleITPC(:,:,t2Tidx,4))/2;
itpcPAAUT(:,:,:,2,2) = (tfSingleITPC(:,:,t2Tidx,1) + tfSingleITPC(:,:,t2Tidx,3))/2;
itpcPAAUT(:,:,:,3,2) = (tfSingleITPC(:,:,t2Tidx,6) + tfSingleITPC(:,:,t2Tidx,8))/2;
itpcPAAUT(:,:,:,4,2) = (tfSingleITPC(:,:,t2Tidx,5) + tfSingleITPC(:,:,t2Tidx,7))/2;

% present vs. absent and attended vs. unattended
for iT = 1:2
    itpcPAT(:,:,:,1,iT) = (itpcPAAUT(:,:,:,1,iT) + itpcPAAUT(:,:,:,2,iT))/2; % present
    itpcPAT(:,:,:,2,iT) = (itpcPAAUT(:,:,:,3,iT) + itpcPAAUT(:,:,:,4,iT))/2; % absent

    itpcAUT(:,:,:,1,iT) = (itpcPAAUT(:,:,:,1,iT) + itpcPAAUT(:,:,:,3,iT))/2; % attended
    itpcAUT(:,:,:,2,iT) = (itpcPAAUT(:,:,:,2,iT) + itpcPAAUT(:,:,:,4,iT))/2; % unattended 
end

% combining across T1 and T2
for iPAAU = 1:4
    itpcPAAU(:,:,:,iPAAU) = (itpcPAAUT(:,:,:,iPAAU,1) + itpcPAAUT(:,:,:,iPAAU,2))/2;
end

itpcPA(:,:,:,1) = (itpcPAAU(:,:,:,1) + itpcPAAU(:,:,:,2))/2; % present
itpcPA(:,:,:,2) = (itpcPAAU(:,:,:,3) + itpcPAAU(:,:,:,4))/2; % absent

itpcAU(:,:,:,1) = (itpcPAAU(:,:,:,1) + itpcPAAU(:,:,:,3))/2; % attended
itpcAU(:,:,:,2) = (itpcPAAU(:,:,:,2) + itpcPAAU(:,:,:,4))/2; % unattended

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
A.stfPADiff(:,:,:,1) = t1SinglePADiff;
A.stfPADiff(:,:,:,2) = t2SinglePADiff;
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

A.itpcAmps = tfSingleITPC;
A.itpcPAAUT = itpcPAAUT;
A.itpcPAT = itpcPAT;
A.itpcAUT = itpcAUT;
A.itpcPAAU = itpcPAAU;
A.itpcPA = itpcPA;
A.itpcAU = itpcAU;

% figures
ytick = 10:10:numel(foi);
xtick = 51:xtickint:numel(toi);
clims = [-0.15 0.15]; 
diffClims = [-0.1 0.1];
hack = plotOrder;
hack(hack>4) = hack(hack>4)+1;

if plotFigs
fH = [];
fH(1) = figure;
set(gcf,'Position',tf9FigPos)
for iTrig = 1:nTrigs
    subplot(2,5,hack(iTrig))
    imagesc(squeeze(nanmean(tfSingleAmps(:,:,:,iTrig),1)),clims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    colormap(parula)
    if iTrig==nTrigs
        xlabel('time (s)')
        ylabel('frequency (Hz)')
    end
    title(trigNames{iTrig})
end
rd_supertitle2('amplitude, all channels mean');

fH(2) = figure;
set(gcf,'Position',tf3FigPos)
attNames = {'attT1','attT2'};
for iAtt = 1:size(tfSingleAmpsAtt,4)
    subplot(1,3,iAtt)
    imagesc(squeeze(nanmean(tfSingleAmpsAtt(:,:,:,iAtt),1)),clims)
    title(attNames{iAtt})
    colormap(parula)
    freezeColors
end
subplot(1,3,3)
imagesc(squeeze(nanmean((tfSingleAmpsAtt(:,:,:,2)-tfSingleAmpsAtt(:,:,:,1)),1)),diffClims)
title('attT2 - attT1')
aH = findall(gcf,'type','axes');
for iAx = 1:numel(aH)
    axes(aH(iAx));
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    xlabel('time (s)')
    ylabel('frequency (Hz)')
end
colormap(cmap)
rd_supertitle2('amplitude, all channels mean');

fH(3) = figure;
set(gcf,'Position',tf9FigPos)
for iPA = 1:size(tfSingleAmpsPA,4)
    subplot(2,4,iPA)
    imagesc(squeeze(nanmean(tfSingleAmpsPA(:,:,:,iPA),1)),clims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    xlabel('time (s)')
    ylabel('frequency (Hz)')
    title(paNames{iPA})
    colormap(parula)
    freezeColors
end
subplot(2,4,5)
imagesc(squeeze(nanmean(t1SinglePADiff,1)),diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (s)')
ylabel('frequency (Hz)')
title(sprintf('T1 %s', paDiffNames))
subplot(2,4,6)
imagesc(squeeze(nanmean(t2SinglePADiff,1)),diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (s)')
ylabel('frequency (Hz)')
title(sprintf('T2 %s', paDiffNames))
subplot(2,4,7)
imagesc(squeeze(nanmean((t2SinglePADiff - t1SinglePADiff),1)),diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (s)')
ylabel('frequency (Hz)')
title(sprintf('T1 vs. T2 %s', paDiffNames))
rd_supertitle2('amplitude, all channels mean');
colormap(cmap)

% 9 squares, attended-unattended
fH(4) = figure;
set(gcf,'Position',tf9SquareFigPos)
% T1/T2 x pres/abs
subplot(3,3,1)
imagesc(squeeze(nanmean((stfPAAUT(:,:,:,1,1)-stfPAAUT(:,:,:,2,1)),1))) % T1-pres-att vs. unatt
ylabel(names{1})
title('T1')
subplot(3,3,2)
imagesc(squeeze(nanmean((stfPAAUT(:,:,:,1,2)-stfPAAUT(:,:,:,2,2)),1))) % T2-pres-att vs. unatt
title('T2')
subplot(3,3,4)
imagesc(squeeze(nanmean((stfPAAUT(:,:,:,3,1)-stfPAAUT(:,:,:,4,1)),1))) % T1-abs-att vs. unatt
ylabel(names{2})
subplot(3,3,5)
imagesc(squeeze(nanmean((stfPAAUT(:,:,:,3,2)-stfPAAUT(:,:,:,4,2)),1))) % T2-abs-att vs. unatt
% ave(T1,T2)
subplot(3,3,3)
imagesc(squeeze(nanmean((stfPAAU(:,:,:,1)-stfPAAU(:,:,:,2)),1))) % pres-att vs. pres-unatt
title('ave(T1,T2)')
subplot(3,3,6)
imagesc(squeeze(nanmean((stfPAAU(:,:,:,3)-stfPAAU(:,:,:,4)),1))) % abs-att vs. abs-unatt
% ave(P,A)
subplot(3,3,7)
imagesc(squeeze(nanmean((stfAUT(:,:,:,1,1)-stfAUT(:,:,:,2,1)),1))) % T1-att vs. T1-unatt 
ylabel(sprintf('ave(%s,\n%s)',names{1},names{2}))
subplot(3,3,8)
imagesc(squeeze(nanmean((stfAUT(:,:,:,1,2)-stfAUT(:,:,:,2,2)),1))) % T2-att vs. T2-unatt 
% ave(all)
subplot(3,3,9)
imagesc(squeeze(nanmean((stfAU(:,:,:,1)-stfAU(:,:,:,2)),1))) % att vs. unatt
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
rd_supertitle2('amplitude, attended vs. unattended')

% 9 squares, present-absent
fH(5) = figure;
set(gcf,'Position',tf9SquareFigPos)
% T1/T2 x att/unatt
subplot(3,3,1)
imagesc(squeeze(nanmean((stfPAAUT(:,:,:,1,1)-stfPAAUT(:,:,:,3,1)),1))) % T1-pres-att vs. abs-att
ylabel('attended')
title('T1')
subplot(3,3,2)
imagesc(squeeze(nanmean((stfPAAUT(:,:,:,1,2)-stfPAAUT(:,:,:,3,2)),1))) % T2-pres-att vs. abs-att
title('T2')
subplot(3,3,4)
imagesc(squeeze(nanmean((stfPAAUT(:,:,:,2,1)-stfPAAUT(:,:,:,4,1)),1))) % T1-pres-unatt vs. abs-unatt
ylabel('unattended')
subplot(3,3,5)
imagesc(squeeze(nanmean((stfPAAUT(:,:,:,2,2)-stfPAAUT(:,:,:,4,2)),1))) % T2-pres-unatt vs. abs-unatt
% ave(T1,T2)
subplot(3,3,3)
imagesc(squeeze(nanmean((stfPAAU(:,:,:,1)-stfPAAU(:,:,:,3)),1))) % pres-att vs. abs-att
title('ave(T1,T2)')
subplot(3,3,6)
imagesc(squeeze(nanmean((stfPAAU(:,:,:,2)-stfPAAU(:,:,:,4)),1))) % pres-unatt vs. abs-unatt
% ave(A,U)
subplot(3,3,7)
imagesc(squeeze(nanmean((stfPAT(:,:,:,1,1)-stfPAT(:,:,:,2,1)),1))) % T1-pres vs. T1-abs 
ylabel('ave(A,U)')
subplot(3,3,8)
imagesc(squeeze(nanmean((stfPAT(:,:,:,1,2)-stfPAT(:,:,:,2,2)),1))) % T2-pres vs. T2-abs 
% ave(all)
subplot(3,3,9)
imagesc(squeeze(nanmean((stfPA(:,:,:,1)-stfPA(:,:,:,2)),1))) % pres vs. abs
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
rd_supertitle2(sprintf('amplitude, %s vs. %s',names{1},names{2}))


if saveFigs
    figPrefix = 'im_wholebrain';
    rd_saveAllFigs(fH, {'timeFreqSingleByCond','timeFreqSingleAtt','timeFreqSinglePA','timeFreqSingleAUDiff','timeFreqSinglePADiff'}, figPrefix, figDir)
end

% itpc
clims = [0 .25];
fH = [];
fH(1) = figure;
set(gcf,'Position',tf9FigPos)
for iTrig = 1:nTrigs
    subplot(2,5,hack(iTrig))
    imagesc(squeeze(nanmean(tfSingleITPC(:,:,:,iTrig),1)),clims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    colormap(parula)
    if iTrig==nTrigs
        xlabel('time (s)')
        ylabel('frequency (Hz)')
    end
    title(trigNames{iTrig})
end
rd_supertitle2('ITPC, all channels mean');

% 9 squares, attended-unattended
fH(2) = figure;
set(gcf,'Position',tf9SquareFigPos)
% T1/T2 x pres/abs
subplot(3,3,1)
imagesc(squeeze(nanmean((itpcPAAUT(:,:,:,1,1)-itpcPAAUT(:,:,:,2,1)),1))) % T1-pres-att vs. unatt
ylabel(names{1})
title('T1')
subplot(3,3,2)
imagesc(squeeze(nanmean((itpcPAAUT(:,:,:,1,2)-itpcPAAUT(:,:,:,2,2)),1))) % T2-pres-att vs. unatt
title('T2')
subplot(3,3,4)
imagesc(squeeze(nanmean((itpcPAAUT(:,:,:,3,1)-itpcPAAUT(:,:,:,4,1)),1))) % T1-abs-att vs. unatt
ylabel(names{2})
subplot(3,3,5)
imagesc(squeeze(nanmean((itpcPAAUT(:,:,:,3,2)-itpcPAAUT(:,:,:,4,2)),1))) % T2-abs-att vs. unatt
% ave(T1,T2)
subplot(3,3,3)
imagesc(squeeze(nanmean((itpcPAAU(:,:,:,1)-itpcPAAU(:,:,:,2)),1))) % pres-att vs. pres-unatt
title('ave(T1,T2)')
subplot(3,3,6)
imagesc(squeeze(nanmean((itpcPAAU(:,:,:,3)-itpcPAAU(:,:,:,4)),1))) % abs-att vs. abs-unatt
% ave(P,A)
subplot(3,3,7)
imagesc(squeeze(nanmean((itpcAUT(:,:,:,1,1)-itpcAUT(:,:,:,2,1)),1))) % T1-att vs. T1-unatt 
ylabel(sprintf('ave(%s,\n%s)',names{1},names{2}))
subplot(3,3,8)
imagesc(squeeze(nanmean((itpcAUT(:,:,:,1,2)-itpcAUT(:,:,:,2,2)),1))) % T2-att vs. T2-unatt 
% ave(all)
subplot(3,3,9)
imagesc(squeeze(nanmean((itpcAU(:,:,:,1)-itpcAU(:,:,:,2)),1))) % att vs. unatt
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
rd_supertitle2('ITPC, attended vs. unattended')

if saveFigs
    figPrefix = 'im_wholebrain';
    rd_saveAllFigs(fH, {'itpcByCond','itpcAUDiff'}, figPrefix, figDir)
end

% multiplots
% channel plot setup
cfg = [];
cfg.layout = layout;
cfg.colormap = cmap;

TFdata.label = data_hdr.label(1:nChannels);
TFdata.dimord = 'chan_freq_time';
TFdata.freq = foi;
TFdata.time = A.stftwinvals;

% A-U amp
vals = stfAU(:,:,:,1)-stfAU(:,:,:,2); % chan x freq x time
TFdata.powspctrm = vals;
cfg.zlim = diffClims*2;

fH = [];
figPos = [1 1 1250 930];
fH(1) = figure('Position',figPos);
ft_multiplotTFR(cfg, TFdata);
rd_supertitle2('amplitude, A-U')

% A-U itpc
vals = itpcAU(:,:,:,1)-itpcAU(:,:,:,2); % chan x freq x time
TFdata.powspctrm = vals;
cfg.zlim = diffClims*2;

figPos = [1 1 1250 930];
fH(2) = figure('Position',figPos);
ft_multiplotTFR(cfg, TFdata);
rd_supertitle2('itpc, A-U')

if saveFigs
    figPrefix = 'immap_wholebrain';
    rd_saveAllFigs(fH, {'timeFreqSingleAUDiff','itpcAUDiff'}, figPrefix, figDir)
end

end

%% save analysis
if saveAnalysis
    save(analysisFileName, 'A') % '-v7.3'
end
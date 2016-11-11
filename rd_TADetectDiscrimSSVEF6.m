function rd_TADetectDiscrimSSVEF6(exptDir, sessionDir, fileBase, analStr, ssvefFreq, trialSelection, respTargetSelection)

% whole brain single trial analysis part 2: fft for select windows and
% frequencies (FOI)

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
        analysisFileName = sprintf('%s/analysis_singleTrials_%s_%s_%sTrials%s_FOI.mat', matDir, fileBase, channelSelectionStr, trialSelection, respTargetSelection);
    otherwise
        savename = sprintf('%s/%s_%s_ssvef_workspace.mat', matDir, fileBase, analStr);
        analysisFileName = sprintf('%s/analysis_singleTrials_%s_%s_%s_%sTrials%s_FOI.mat', matDir, fileBase, analStr, channelSelectionStr, trialSelection, respTargetSelection);
end

%% Get the data
load(savename)

%% Update behav
behav = behavior(behav);

%% Settings after loading the data
saveAnalysis = 1;
saveFigs = 0;
plotFigs = 0;

excludeTrialsFt = 1;
excludeSaturatedEpochs = 0;

channels = 1:157;

%% Plotting setup
if plotFigs
load parula
cmap = flipud(lbmap(64,'RedBlue'));
attNames = {'att','unatt'};
set(0,'defaultLineLineWidth',1)

% load data header for plotting topologies
load data/data_hdr.mat
cfg = [];
layout = ft_prepare_layout(cfg, data_hdr);
end

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
            analysisFileName = sprintf('%s/analysis_singleTrials_%s_ft_%s_%sTrials%s_FOI.mat', matDir, fileBase, channelSelectionStr, trialSelection, respTargetSelection);
        otherwise
            analysisFileName = sprintf('%s/analysis_singleTrials_%s_%s_ft_%s_%sTrials%s_FOI.mat', matDir, fileBase, analStr, channelSelectionStr, trialSelection, respTargetSelection);
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

A.nTrialsCond = nTrialsCond;

%% 600 ms pre-stimulus
%% Select time windows of interest
twin = [-600 0];
t1Tidx = find(t==eventTimes(3)+twin(1)):find(t==eventTimes(3)+twin(2));
t2Tidx = find(t==eventTimes(4)+twin(1)):find(t==eventTimes(4)+twin(2));

preAU = [];
% window before T1, separated by attention condition
valsAtt = trigMean(t1Tidx,:,:,[1 3 5 7]); % att T1
valsUnatt = trigMean(t1Tidx,:,:,[2 4 6 8]); % att T2
preAU{1}(:,:,:,1) = reshape(valsAtt,numel(t1Tidx),nChannels,nTrialsPerCond*4);
preAU{1}(:,:,:,2) = reshape(valsUnatt,numel(t1Tidx),nChannels,nTrialsPerCond*4);

% window before T2, T1 absent only, separated by attention condition
valsAtt = trigMean(t2Tidx,:,:,[4 8]); % att T2
valsUnatt = trigMean(t2Tidx,:,:,[3 7]); % att T1
preAU{2}(:,:,:,1) = reshape(valsAtt,numel(t1Tidx),nChannels,nTrialsPerCond*2);
preAU{2}(:,:,:,2) = reshape(valsUnatt,numel(t1Tidx),nChannels,nTrialsPerCond*2);

A.preAU600 = preAU;

%% FFT
ampsAll = [];
preAUTAmps = [];
for iT = 1:2
    % taper
    twinvals = twin(1):twin(end);
    nSamples = length(twinvals);
    vals = preAU{iT};
    sz = size(vals);
    taper = hanning(nSamples);
    vals = vals.*repmat(taper,1,sz(2),sz(3),sz(4));
    
    % fft
    nfft = 1000; %2^nextpow2(nSamples); % Next power of 2 from length of y
    Y = fft(vals,nfft)/nSamples; % Scale by number of samples
    f = Fs/2*linspace(0,1,nfft/2+1); % Fs/2 is the maximum frequency that can be measured
    amps = 2*abs(Y(1:nfft/2+1,:,:,:)); % Multiply by 2 since only half the energy is in the positive half of the spectrum?
    
    ampsAll{iT} = amps;
    
    % average across trials
    preAUTAmps(:,:,:,iT) = squeeze(nanmean(amps,3));
end

% calculate z-score across trials
ampsAllC = cat(3,ampsAll{1},ampsAll{2}); % combine T1 and T2
sz = size(ampsAllC);
ampsAllC = reshape(ampsAllC,sz(1),sz(2),sz(3)*sz(4));
m = nanmean(ampsAllC,3);
sd = nanstd(ampsAllC,0,3);

preAUTAmpsZ = [];
for iT = 1:2
    ampsAllZ = (ampsAll{iT}-repmat(m,1,1,size(ampsAll{iT},3),size(ampsAll{iT},4)))./...
        repmat(sd,1,1,size(ampsAll{iT},3),size(ampsAll{iT},4));
    preAUTAmpsZ(:,:,:,iT) = squeeze(nanmean(ampsAllZ,3));
end

A.preAUTAmps600 = preAUTAmps;
A.preAUTAmpsZ600 = preAUTAmpsZ;

%% Alpha
freqRange = [8 12];
fIdx = f>=freqRange(1) & f<=freqRange(2);
preAUT = squeeze(mean(preAUTAmps(fIdx,:,:,:),1));
preAUTZ = squeeze(mean(preAUTAmpsZ(fIdx,:,:,:),1));

if plotFigs
clims = [10 70];
diffClims = [-5 5];
fH = [];
fH(1) = figure('Position',[360 280 450 630]);
for iT = 1:2
    for iAU = 1:2
        subplot(3,2,iT+2*(iAU-1))
        str = sprintf('T%d %s',iT, attNames{iAU});
        vals = preAUT(:,iAU,iT)';
        ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
        set(gca,'CLim',clims)
        colormap(parula)
        freezeColors
    end
    subplot(3,2,iT+4)
    str = sprintf('T%d att-unatt',iT);
    vals = preAUT(:,1,iT)'-preAUT(:,2,iT)';
    ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
    set(gca,'CLim',diffClims)
    colormap(cmap)
end
rd_supertitle2(sprintf('f = [%d %d] Hz',freqRange(1),freqRange(2)))
end

A.alpha.twin = twin;
A.alpha.freqRange = freqRange;
A.alpha.preAUT = preAUT;
A.alpha.preAUTZ = preAUTZ;

%% 200 ms pre-stimulus
%% Select time windows of interest
twin = [-200 0];
t1Tidx = find(t==eventTimes(3)+twin(1)):find(t==eventTimes(3)+twin(2));
t2Tidx = find(t==eventTimes(4)+twin(1)):find(t==eventTimes(4)+twin(2));

preAU = [];
% window before T1, separated by attention condition
valsAtt = trigMean(t1Tidx,:,:,[1 3 5 7]); % att T1
valsUnatt = trigMean(t1Tidx,:,:,[2 4 6 8]); % att T2
preAU{1}(:,:,:,1) = reshape(valsAtt,numel(t1Tidx),nChannels,nTrialsPerCond*4);
preAU{1}(:,:,:,2) = reshape(valsUnatt,numel(t1Tidx),nChannels,nTrialsPerCond*4);

% window before T2, T1 absent only, separated by attention condition
valsAtt = trigMean(t2Tidx,:,:,[4 8]); % att T2
valsUnatt = trigMean(t2Tidx,:,:,[3 7]); % att T1
preAU{2}(:,:,:,1) = reshape(valsAtt,numel(t1Tidx),nChannels,nTrialsPerCond*2);
preAU{2}(:,:,:,2) = reshape(valsUnatt,numel(t1Tidx),nChannels,nTrialsPerCond*2);

A.preAU200 = preAU;

%% FFT
preAUTAmps = [];
for iT = 1:2
    % taper
    twinvals = twin(1):twin(end);
    nSamples = length(twinvals);
    vals = preAU{iT};
    sz = size(vals);
    taper = hanning(nSamples);
    vals = vals.*repmat(taper,1,sz(2),sz(3),sz(4));
    
    % fft
    nfft = 1000; % Next power of 2 from length of y
    Y = fft(vals,nfft)/nSamples; % Scale by number of samples
    f = Fs/2*linspace(0,1,nfft/2+1); % Fs/2 is the maximum frequency that can be measured
    amps = 2*abs(Y(1:nfft/2+1,:,:,:)); % Multiply by 2 since only half the energy is in the positive half of the spectrum?
    
    ampsAll{iT} = amps;
    
    % average across trials
    preAUTAmps(:,:,:,iT) = squeeze(nanmean(amps,3));
end

% calculate z-score across trials
ampsAllC = cat(3,ampsAll{1},ampsAll{2}); % combine T1 and T2
sz = size(ampsAllC);
ampsAllC = reshape(ampsAllC,sz(1),sz(2),sz(3)*sz(4));
m = nanmean(ampsAllC,3);
sd = nanstd(ampsAllC,0,3);

preAUTAmpsZ = [];
for iT = 1:2
    ampsAllZ = (ampsAll{iT}-repmat(m,1,1,size(ampsAll{iT},3),size(ampsAll{iT},4)))./...
        repmat(sd,1,1,size(ampsAll{iT},3),size(ampsAll{iT},4));
    preAUTAmpsZ(:,:,:,iT) = squeeze(nanmean(ampsAllZ,3));
end

A.preAUTAmps200 = preAUTAmps;
A.preAUTAmpsZ200 = preAUTAmpsZ;

%% SSVEF 30
freqRange = [30 30];
fIdx = f>=freqRange(1) & f<=freqRange(2);
preAUT = squeeze(mean(preAUTAmps(fIdx,:,:,:),1));
preAUTZ = squeeze(mean(preAUTAmpsZ(fIdx,:,:,:),1));

if plotFigs
clims = [0 40];
diffClims = [-4 4];
fH(2) = figure('Position',[360 280 450 630]);
for iT = 1:2
    for iAU = 1:2
        subplot(3,2,iT+2*(iAU-1))
        str = sprintf('T%d %s',iT, attNames{iAU});
        vals = preAUT(:,iAU,iT)';
        ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
        set(gca,'CLim',clims)
        colormap(parula)
        freezeColors
    end
    subplot(3,2,iT+4)
    str = sprintf('T%d att-unatt',iT);
    vals = preAUT(:,1,iT)'-preAUT(:,2,iT)';
    ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
    set(gca,'CLim',diffClims)
    colormap(cmap)
end
rd_supertitle2(sprintf('f = [%d %d] Hz',freqRange(1),freqRange(2)))
end

A.ssvef30.twin = twin;
A.ssvef30.freqRange = freqRange;
A.ssvef30.preAUT = preAUT;
A.ssvef30.preAUTZ = preAUTZ;

%% SSVEF 40
freqRange = [40 40];
fIdx = f>=freqRange(1) & f<=freqRange(2);
preAUT = squeeze(mean(preAUTAmps(fIdx,:,:,:),1));
preAUTZ = squeeze(mean(preAUTAmpsZ(fIdx,:,:,:),1));

if plotFigs
clims = [0 40];
diffClims = [-4 4];
fH(3) = figure('Position',[360 280 450 630]);
for iT = 1:2
    for iAU = 1:2
        subplot(3,2,iT+2*(iAU-1))
        str = sprintf('T%d %s',iT, attNames{iAU});
        vals = preAUT(:,iAU,iT)';
        ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
        set(gca,'CLim',clims)
        colormap(parula)
        freezeColors
    end
    subplot(3,2,iT+4)
    str = sprintf('T%d att-unatt',iT);
    vals = preAUT(:,1,iT)'-preAUT(:,2,iT)';
    ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
    set(gca,'CLim',diffClims)
    colormap(cmap)
end
rd_supertitle2(sprintf('f = [%d %d] Hz',freqRange(1),freqRange(2)))
end

A.ssvef40.twin = twin;
A.ssvef40.freqRange = freqRange;
A.ssvef40.preAUT = preAUT;
A.ssvef40.preAUTZ = preAUTZ;

%% Broadband
freqRange = [70 100];
fIdx = f>=freqRange(1) & f<=freqRange(2);
preAUT = squeeze(mean(preAUTAmps(fIdx,:,:,:),1));
preAUTZ = squeeze(mean(preAUTAmpsZ(fIdx,:,:,:),1));

if plotFigs
clims = [0 20];
diffClims = [-2 2];
fH(4) = figure('Position',[360 280 450 630]);
for iT = 1:2
    for iAU = 1:2
        subplot(3,2,iT+2*(iAU-1))
        str = sprintf('T%d %s',iT, attNames{iAU});
        vals = preAUT(:,iAU,iT)';
        ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
        set(gca,'CLim',clims)
        colormap(parula)
        freezeColors
    end
    subplot(3,2,iT+4)
    str = sprintf('T%d att-unatt',iT);
    vals = preAUT(:,1,iT)'-preAUT(:,2,iT)';
    ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
    set(gca,'CLim',diffClims)
    colormap(cmap)
end
rd_supertitle2(sprintf('f = [%d %d] Hz',freqRange(1),freqRange(2)))
end

A.broadband.twin = twin;
A.broadband.freqRange = freqRange;
A.broadband.preAUT = preAUT;
A.broadband.preAUTZ = preAUTZ;

%% Save figs
if saveFigs
    figPrefix = 'map_wholebrain';
    rd_saveAllFigs(fH, {'pre600_8-12Hz','pre200_30Hz','pre200_40Hz','pre200_70-100Hz'}, figPrefix, figDir)
end

%% save analysis
if saveAnalysis
    save(analysisFileName, 'A') 
end

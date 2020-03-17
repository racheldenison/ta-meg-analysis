function corrTS = rd_timeSegmentReliability(sessionDir)

%% i/o
exptDir = pathToTANoise('MEG');
% sessionDir = 'R0817_20171212';

dataFile = dir(sprintf('%s/%s/mat/*condData.mat', exptDir, sessionDir));
dataFileName = sprintf('%s/%s/mat/%s', exptDir, sessionDir, dataFile.name);
channelsFile = dir(sprintf('%s/%s/mat/channels*.mat', exptDir, sessionDir));
channelsFileName = sprintf('%s/%s/mat/%s', exptDir, sessionDir, channelsFile.name);
% figDir = sprintf('%s/%s/figures/ebi_ft', exptDir, sessionDir);

nCh = 5;

% saveFigs = 0;

%% load data
fprintf('\n\n%s', sessionDir)
fprintf('\nLoading data ... ')
load(dataFileName)
fprintf('done')

% load channels
C = load(channelsFileName);
channels = C.channelsRanked(1:nCh);

% load data header for plotting topologies
load data/data_hdr.mat

%% setup
t = D.t;
Fs = D.Fs;
eventTimes = D.eventTimes;
data = D.condData;

nT = numel(t);
sz = size(data);
nCue = sz(1);
nT1 = sz(2);
nT2 = sz(3);

% cueNames = {'precue T1','precue T2','neutral'};
cueNames = {'precue T1','precue T2'};

segmentLength = 100;
nPerm = 100;

%% format data
if ~iscell(data)
    sz = size(data);
    nCue = sz(4);
    nT1 = sz(5);
    nT2 = sz(6);
    data0 = data;
    data = [];
    for iCue = 1:nCue
        for iT1 = 1:nT1
            for iT2 = 1:nT2
                data{iCue,iT1,iT2} = data0(:,:,:,iCue,iT1,iT2);
            end
        end
    end
end

%% organize data by cue
cueData = [];
for iCue = 1:nCue
    cueData{iCue} = [];
    for iT1 = 1:nT1
        for iT2 = 1:nT2
            cueData{iCue} = cat(3, cueData{iCue}, data{iCue,iT1,iT2}(:,channels,:));
        end
    end
end

nTrialsPerCue = size(cueData{1},3);

%% split half time segment correlation
corrTS = [];
for iCue = 1:nCue
    fprintf('\n%s\n', cueNames{iCue})
    for iP = 1:nPerm
        fprintf('.')
        if mod(iP,10)==0, fprintf('%d',iP), end
        idx = mod(randperm(nTrialsPerCue),2)==1;
        
        vals1 = nanmean(cueData{iCue}(:,:,idx),3);
        vals2 = nanmean(cueData{iCue}(:,:,~idx),3);
        
        corrvals = nan(nT,nCh);
        for iT = segmentLength/2+1:nT-segmentLength/2
            tidx = iT-segmentLength/2:iT+segmentLength/2;
            corrvals(iT,:) = diag(corr(vals1(tidx,:),vals2(tidx,:)));
        end
        
        corrTS{iCue}(:,:,iP) = corrvals;
    end
end

figure
hold on
for iCue = 1:nCue
    plot(t, mean(mean(corrTS{iCue},3),2))
end
title(sprintf('%s', und2space(sessionDir)))


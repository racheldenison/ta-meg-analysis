function distTS = rd_timeSegmentReliabilityDist(sessionDir)

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
% load data/data_hdr.mat

%% setup
t = D.t;
data = D.condData;

nT = numel(t);
sz = size(data);
nCue = sz(1);
nT1 = sz(2);
nT2 = sz(3);

% cueNames = {'precue T1','precue T2','neutral'};
cueNames = {'precue T1','precue T2'};

segmentLength = 50;

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

%% time segment distance
% each trial segment is a point, quantify distance among points
% pdist treats row vectors as obesrvations
distTS = [];
for iCue = 1:nCue
    fprintf('\n%s\n', cueNames{iCue})
    
    distvals = nan(nT,nCh);
    for iCh = 1:nCh
        fprintf('.')
        for iT = segmentLength/2+1:nT-segmentLength/2
            tidx = iT-segmentLength/2:iT+segmentLength/2;
            vals = squeeze(cueData{iCue}(tidx,iCh,:));
            valsz = zscore(vals); % zscore to emphasize pattern across time
            nanidx = any(isnan(vals));
            dz = pdist(valsz(:,~nanidx)'); % Euclidean distances between pairs of trials (use squareform to plot)

            distvals(iT,iCh) = mean(dz(dz~=0));
        end
    end
    
    distTS{iCue} = distvals;
end

figure
hold on
for iCue = 1:nCue
    plot(t, mean(distTS{iCue},2))
end
title(sprintf('%s', und2space(sessionDir)))


% rd_segmentSqd.m
%
% Chunk data into smaller segments. Segments data between runs, based on 
% triggers.
%
% Rachel Denison
% September 2014

%% setup
exptDir = '/Local/Users/denison/Data/TAPilot/MEG';
sessionDir = 'R0890_20140806';
dataFile = 'R0890_TAPilot_8.06.14.sqd';
fileName = sprintf('%s/%s/%s', exptDir, sessionDir, dataFile);

trigChan = 166; % blank blocks

nTrigsPerRun = 9;
nRuns = 15;
nRunsPerSegment = 3;
trialDur = 7;

%% display sqd info
info = sqdread(fileName,'info')

%% read in all triggers
triggers = all_trigger(fileName, trigChan);
nTrigs = size(triggers,1);

%% check and inspect triggers
nTrigsExpected = nTrigsPerRun*nRuns;

if nTrigs~=nTrigsExpected
    fprintf('\nFound %d triggers, but expected %d!\n\n', nTrigs, nTrigsExpected)
    
    figure
    plot(triggers(:,1),triggers(:,2),'.')
    title('original triggers')
end

%% exclude stray triggers
% manual (different for each subject)
if nTrigs>nTrigsExpected
    excludedTrigIdxs = nTrigsPerRun*6 + 1;
    triggers(excludedTrigIdxs,:) = [];
    nTrigs = size(triggers,1);

    fprintf('\nNow there are %d triggers (%d expected)\n\n', nTrigs, nTrigsExpected)
    figure
    plot(triggers(:,1),triggers(:,2),'.')
    title('after trigger exclusion')
end

%% specify data segments
nTrigsPerSegment = nRunsPerSegment*nTrigsPerRun;
segmentFirstTrialIdxs = 1:nTrigsPerSegment:nTrigs;
segmentLastTrialIdxs = [segmentFirstTrialIdxs(2:end)-1 nTrigs];
segmentFirstTrialStartTimes = triggers(segmentFirstTrialIdxs,1);
segmentLastTrialEndTimes = triggers(segmentLastTrialIdxs,1) + trialDur;

% % divide segments halfway through last trial end and first trial start of
% % subsequent segments
segmentDividePoints = mean([segmentLastTrialEndTimes(1:end-1) segmentFirstTrialStartTimes(2:end)],2)';
segmentStartTimes = [int64(1) segmentDividePoints];
segmentEndTimes = [segmentDividePoints+1 int64(info.SamplesAvailable)];
nSegments = numel(segmentStartTimes);

%% read data segments and write new sqd file
for iSegment = 1:nSegments
	segment = sqdread(fileName, 'Samples', [segmentStartTimes(iSegment) segmentEndTimes(iSegment)]);
    
    segmentFileName = sprintf('%s_segment%d.sqd', fileName(1:end-4), iSegment);
    sqdwrite(fileName, segmentFileName, 'Data', segment);
end




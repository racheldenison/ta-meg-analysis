function nRuns = rd_segmentSqd(fileName)
%
% function nRuns = rd_segmentSqd(fileName)
%
% Chunk data into smaller segments. Segments data between runs, based on 
% triggers.
%
% Rachel Denison
% September 2014

%% setup
% file in
% exptDir = '/Local/Users/denison/Data/TAPilot/MEG';
% sessionDir = 'R0890_20140806';
% dataFile = 'R0890_TAPilot_8.06.14.sqd';
% fileName = sprintf('%s/%s/%s', exptDir, sessionDir, dataFile);

% file out
segmentLabel = 'run';

% trigger
trigChan = 166; % 166=blank blocks % for dividing runs

% segmentation options
segmentationOption = 'selectData'; % 'splitData' or 'selectData'
segmentCushion = 5; % used only if 'selectData'

% experiment info
nTrigsPerRun = 9; % TAPilot: 9
nRuns = 9; % TADetectDiscrim: 14; TAPilot: 18
nRunsPerSegment = 1;
segmentOutNums = []; % [91 1:14]; % default is [], numbers 1:nRuns
trialDur = 7;

%% display sqd info
info = sqdread(fileName,'info')

%% read in all triggers
triggers = all_trigger(fileName, trigChan);
nTrigs = size(triggers,1);

%% check and inspect triggers
nTrigsExpected = nTrigsPerRun*nRuns;

if nTrigs~=nTrigsExpected
    fprintf('Found %d triggers, but expected %d!\n', nTrigs, nTrigsExpected)
    
    figure
    plot(triggers(:,1),triggers(:,2),'.')
    title('original triggers')
    
    fprintf('\nExiting to allow manual adjustments ...\n')
    return
end

%% exclude stray triggers
% manual (different for each subject)
if nTrigs>nTrigsExpected
%     excludedTrigIdxs = nTrigsPerRun*6 + 1;
    excludedTrigIdxs = [1:8]; % TADeDi R1021_20151120 r1-5: [1 20:26]; r6-14: [1:5]; R1026_20151211: [82:84 112:115]; R0582_20151211 r1: [1:8], r2-14: [10:11 21 94:97]
    triggers(excludedTrigIdxs,:) = [];
    nTrigs = size(triggers,1);

    fprintf('\nNow there are %d triggers (%d expected)\n\n', nTrigs, nTrigsExpected)
    figure
    plot(triggers(:,1),triggers(:,2),'.')
    title('after trigger exclusion')
end

%% specify data segments
Fs = info.SampleRate;

nTrigsPerSegment = nRunsPerSegment*nTrigsPerRun;
segmentFirstTrialIdxs = 1:nTrigsPerSegment:nTrigs;
segmentLastTrialIdxs = [segmentFirstTrialIdxs(2:end)-1 nTrigs];
segmentFirstTrialStartTimes = triggers(segmentFirstTrialIdxs,1);
segmentLastTrialEndTimes = triggers(segmentLastTrialIdxs,1) + trialDur*Fs;

switch segmentationOption
    case 'splitData'
        % % divide segments halfway through last trial end and first trial start of
        % % subsequent segments
        segmentDividePoints = mean([segmentLastTrialEndTimes(1:end-1) segmentFirstTrialStartTimes(2:end)],2)';
        segmentStartTimes = [int64(1) segmentDividePoints];
        segmentEndTimes = [segmentDividePoints+1 int64(info.SamplesAvailable)];
    case 'selectData'
        segmentStartTimes = segmentFirstTrialStartTimes - segmentCushion*Fs;
        segmentEndTimes = segmentLastTrialEndTimes + segmentCushion*Fs;
    otherwise
        error('segmentationOption not found')
end
nSegments = numel(segmentStartTimes);
if isempty(segmentOutNums)
    segmentOutNums = 1:nSegments;
end

%% read data segments and write new sqd file
for iSegment = 1:nSegments
    % read segment data from original file
	segment = sqdread(fileName, 'Samples', [segmentStartTimes(iSegment) segmentEndTimes(iSegment)]);
    
    % new segment file name
    segmentOutNum = segmentOutNums(iSegment);
    segmentFileName = sprintf('%s_%s%02d.sqd', fileName(1:end-4), segmentLabel, segmentOutNum);
    
    % write segment file 
    fprintf('Writing sqd: %s %d\n', segmentLabel, segmentOutNum) 
    sqdwrite(fileName, segmentFileName, 'Data', segment);
    
    % check triggers in new file
    rd_checkTriggers(segmentFileName,[],0);
    title(sprintf('%s %d', segmentLabel, segmentOutNum)) 
end





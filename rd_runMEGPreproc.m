% rd_runMEGPreproc.m

%% setup
% exptDir = '/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG';
% exptDir = '/Local/Users/denison/Data/TANoise/MEG';
% exptDir = '/Local/Users/denison/Data/TA2/MEG';
exptDir = pathToTA2('MEG');

sessionDir = 'R1507_20190621';
fileBase = 'R1507_TA2_6.21.19';

renameFiles = false;
runsToRename = 1:12;

segmentDataFile = false;

dataDir = sprintf('%s/%s', exptDir, sessionDir);
preprocDir = sprintf('%s/preproc', dataDir);
figDir = sprintf('%s/%s/%s', preprocDir, 'figures');

inspectData = false;

%% rename files if needed
if renameFiles
    renameRunFiles(sessionDir, runsToRename) 
end

%% make the preproc dir if it doesn't exist
if ~exist(preprocDir,'dir')
    mkdir(preprocDir)
end

%% segment only if needed
runFiles = dir(sprintf('%s/%s*.sqd', preprocDir, fileBase));
if isempty(runFiles)
    if segmentDataFile
        %% segment original sqd into runs
        dataFile = sprintf('%s/%s.sqd', dataDir, fileBase);
        
        % check settings in rd_segmentSqd before running!
        nRuns = rd_segmentSqd(dataFile);
        runs = 1:nRuns
    end
    
    %% move run files into preproc directory
    runFiles = dir(sprintf('%s/*run*.sqd', dataDir));
    nRuns = numel(runFiles);
    runs = 1:nRuns
    
    for iRun = 1:nRuns
        movefile(sprintf('%s/%s', dataDir, runFiles(iRun).name), preprocDir)
    end
else
    % we have done preprocessing before, so find the number of runs
%     runTag = rd_getTag(runFiles(end).name,'run');
%     nRuns = str2num(runTag(3:4));
    nRuns = numel(runFiles);
    runs = 1:nRuns
end

%% view data
if inspectData
    % % just from run 1
    run1Data = sqdread(sprintf('%s/%s', preprocDir, runFiles(1).name));
    run1Data  = run1Data(:,1:157)';
    
    srate = 1000;
    windowSize = [1 5 2560 1392];
    eegplot(run1Data,'srate',srate,'winlength',20,'dispchans',80,'position',windowSize);
end

%% manually set bad channels
badChannels = []; % in matlab 1-indexing

%% run preproc for each run
for iRun = 1:nRuns
    run = runs(iRun);
%     runFile = sprintf('%s/%s_run%02d.sqd', preprocDir, fileBase, run);
    runFile = sprintf('%s/%s', preprocDir, runFiles(iRun).name);
    preprocFileName = rd_MEGPreproc(runFile, figDir, badChannels);
end

%% combine run files into preprocessed sqd
analStr = rd_getTag(preprocFileName);
outFileName = sprintf('%s_%s.sqd', fileBase, analStr);
outfile = rd_combineSqd(preprocDir, outFileName, analStr)

%% view triggers for the combined file
rd_checkTriggers(outfile);


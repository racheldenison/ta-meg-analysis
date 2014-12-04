% rd_runMEGPreproc.m

%% setup
exptDir = '/Local/Users/denison/Data/TAPilot/MEG';
sessionDir = 'R0890_20140806';
fileBase = 'R0890_TAPilot_8.06.14';

dataDir = sprintf('%s/%s', exptDir, sessionDir);
preprocDir = sprintf('%s/preproc', dataDir);
figDir = sprintf('%s/%s/%s', preprocDir, 'figures');

%% make the preproc dir if it doesn't exist
if ~exist(preprocDir,'dir')
    mkdir(preprocDir)
end

%% segment only if needed
runFiles = dir(sprintf('%s/%s*.sqd', preprocDir, fileBase));
if isempty(runFiles)
    %% segment original sqd into runs
    dataFile = sprintf('%s/%s.sqd', dataDir, fileBase);
    
    % check settings in rd_segmentSqd before running!
    nRuns = rd_segmentSqd(dataFile);
    runs = 1:nRuns
    
    %% move run files into preproc directory
    runFiles = dir(sprintf('%s/*run*.sqd', dataDir));
    for iRun = 1:nRuns
        movefile(sprintf('%s/%s', dataDir, runFiles(iRun).name), preprocDir)
    end
else
    % we have done preprocessing before, so find the number of runs
    runTag = rd_getTag(runFiles(end).name,'run');
    nRuns = str2num(runTag(3:4));
    runs = 1:nRuns
end

%% run preproc for each run
for iRun = 1:nRuns
    run = runs(iRun);
    runFile = sprintf('%s/%s_run%02d.sqd', preprocDir, fileBase, run);
    preprocFileName = rd_MEGPreproc(runFile, figDir);
end

%% combine run files into preprocessed sqd
analStr = rd_getTag(preprocFileName);
outFileName = sprintf('%s_%s.sqd', fileBase, analStr);
outfile = rd_combineSqd(preprocDir, outFileName, analStr)

%% view triggers for the combined file
rd_checkTriggers(outfile);


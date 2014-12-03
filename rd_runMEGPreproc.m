% rd_runMEGPreproc.m

%% setup
exptDir = '/Local/Users/denison/Data/TAPilot/MEG';
sessionDir = 'R0817_20140820';
fileBase = 'R0817_TAPilot_8.20.14';
analStr = 'e';

dataDir = sprintf('%s/%s', exptDir, sessionDir);
preprocDir = sprintf('%s/preproc', dataDir);
figDir = sprintf('%s/%s', preprocDir, 'figures');
preprocFileName = sprintf('%s_%s.sqd', fileBase, analStr);

%% make the preproc dir if it doesn't exist
if ~exist(preprocDir,'dir')
    mkdir(preprocDir)
end

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

%% run preproc for each run
for iRun = 1:nRuns
    run = runs(iRun);
    runFile = sprintf('%s/%s_run%02d.sqd', preprocDir, fileBase, run);
    rd_MEGPreproc(runFile, figDir);
end

%% combine run files into preprocessed sqd
outfile = rd_combineSqd(preprocDir, preprocFileName, analStr)

%% view triggers for the combined file
rd_checkTriggers(outfile);


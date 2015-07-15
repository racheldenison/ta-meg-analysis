% rd_makeTrialsRejected.m

%% setup
exptDir = '/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG';
sessionDir = 'R0817_20150526';
fileBase = 'R0817_TADeDi_5.26.15';
preprocDir = sprintf('%s/%s/preproc', exptDir, sessionDir);
matDir = sprintf('%s/%s/mat', exptDir, sessionDir);

prepFiles = dir([preprocDir '/*prepCleanData.mat']);

runsPerFile = [5 1 3 3 3];
trialsPerRun = 41;

specialAdjustments = 1;

%% check files
for iFile = 1:numel(prepFiles)
    fprintf('%d runs\t', runsPerFile(iFile))
    disp(prepFiles(iFile).name)
end
ok = input('\nAre these the files in the correct order? [y/n]','s');
if ~strcmp(ok,'y')
    fprintf('\nok, quitting ...\n')
    break
end

%% get rejected trials from each prep file
for iFile = 1:numel(prepFiles)
    prepFile = prepFiles(iFile).name;
    load(sprintf('%s/%s', preprocDir, prepFile))
    tr{iFile} = cleanPrepData.trials_rejected(:,1);
    if iFile==1
        startingTrial(iFile) = 1;
        endingTrial(iFile) = runsPerFile(iFile)*trialsPerRun;
    else
        startingTrial(iFile) = endingTrial(iFile-1)+1;
        endingTrial(iFile) = startingTrial(iFile) + runsPerFile(iFile)*trialsPerRun - 1;
    end
end

%% special adjustments (oy)
if specialAdjustments
    ok = input('\nAre you using the right adjustments? [y/n]','s');
    if ~strcmp(ok,'y')
        fprintf('\nok, quitting ...\n')
        break
    end
    
    % R0817_20150526
    iswrong = tr{5}>44;
    tr{5}(iswrong) = tr{5}(iswrong)-5;
    tr{5} = unique(tr{5});
end

%% construct trials_rejected
trials_rejected = [];
for iFile = 1:numel(prepFiles)
    trials_rejected = [trials_rejected; tr{iFile} + startingTrial(iFile) - 1];
end

%% save trials_rejected
save(sprintf('%s/trials_rejected.mat', matDir), 'trials_rejected')


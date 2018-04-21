function renameRunFiles(sessionName, runs)

%% example inputs
% sessionName = 'R0817_20171213';
% runs = 1:12;

%% setup
% exptDir = '/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG';
exptDir = '/Local/Users/denison/Data/TANoise/MEG';
exptName = 'TANoise';
exptNameJeff = 'TA_Noise';
dataDir = sprintf('%s/%s', exptDir, sessionName);

C = strsplit(sessionName,'_');
subject = C{1};
dateStr = C{2};
dateStrJeff = sprintf('%s.%s.%s', dateStr(5:6), dateStr(7:8), dateStr(3:4));

%% rename
for iRun = 1:numel(runs)
    fileNameJeff = sprintf('%s/%s_%s_Block%d_%s.sqd', dataDir, subject, exptNameJeff, runs(iRun), dateStrJeff);
    fileNameNew = sprintf('%s/%s_%s_%s_run%02d.sqd', dataDir, subject, exptName, dateStrJeff, runs(iRun));
    
    if iRun==1
        fprintf('\nExample old file:\n%s\n', fileNameJeff)
        fprintf('\nExample new file:\n%s\n', fileNameNew)
        answer = input(sprintf('\nOk to rename all %d runs? [y/n] ', numel(runs)),'s');
    end
    if ~strcmp(answer,'y')
        fprintf('\nNot renaming, exiting ...\n')
        break
    end
    
    movefile(fileNameJeff, fileNameNew)
end


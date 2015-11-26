% rd_runTADetectDiscrimSSVEF2.m

%% setup
exptDir = '/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG';
analStr = 'ebi'; % '', 'ebi', etc.
ssvefFreqs = [30 40];
nTopChs = [1 5]; % 1, 5, etc.

% subjects = {'R0817_20150504', 'R0973_20150727', 'R0974_20150728', ...
%     'R0861_20150813', 'R0504_20150805', 'R0983_20150813', ...
%     'R0898_20150828', 'R0436_20150904', 'R0988_20150904'};
subjects = {'R1021_20151120'};
nSubjects = numel(subjects);

%% run analysis
for iSubject = 1:nSubjects
    % get fileBase
    sessionDir = subjects{iSubject};
    fprintf('%s\n',sessionDir)
    
    sqdFile = dir(sprintf('%s/%s/*_%s.sqd', exptDir, sessionDir, analStr));
    if numel(sqdFile)~=1
        error('More or fewer than 1 matching data file\n%s/%s/*_%s.sqd', exptDir, sessionDir, analStr);
    end
    [~,didx] = getTag(sqdFile.name,'_');
    fileBase = sqdFile.name(1:didx-1);
    
    % run freq/top channels combos
    for ssvefFreq = ssvefFreqs
        for nTopChannels = nTopChs
            rd_TADetectDiscrimSSVEF2(exptDir, sessionDir, fileBase, ...
                analStr, ssvefFreq, nTopChannels);
            close all;
        end
    end
end
fprintf('done.\n')
    
    
    
    
    
    
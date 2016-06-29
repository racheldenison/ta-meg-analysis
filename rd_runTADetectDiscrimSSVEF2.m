% rd_runTADetectDiscrimSSVEF2.m

%% setup
exptDir = '/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG';
analStr = 'ebi'; % '', 'ebi', etc.
ssvefFreqs = [30 40];
nTopChs = 5; % 1, 5, [1 5], etc.
iqrThreshs = [];
weightChannels = 0;
trialSelection = 'detectHit'; % 'all','validCorrect','detectHit','detectMiss','detectFA','detectCR','discrimCorrect','discrimIncorrect'

subjects = {'R0817_20150504', 'R0973_20150727', 'R0974_20150728', ...
    'R0861_20150813', 'R0504_20150805', 'R0983_20150813', ...
    'R0898_20150828', 'R0436_20150904', 'R1018_20151118', ...
    'R1019_20151118','R1021_20151120','R1026_20151211', ...
    'R0852_20151211','R1027_20151216','R1028_20151216',...
    'R1029_20151222'}; % N=16
% subjects = subjects(8:end);
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
        
        if ~isempty(nTopChs) && ~isempty(iqrThreshs)
            error('set either nTopChs or iqrThreshs to empty')
        else
            if ~isempty(nTopChs)
                for nTopChannels = nTopChs
                    rd_TADetectDiscrimSSVEF2(exptDir, sessionDir, fileBase, ...
                        analStr, ssvefFreq, nTopChannels, [], weightChannels, trialSelection);
                    close all;
                end
            elseif ~isempty(iqrThreshs)
                for iqrThresh = iqrThreshs
                    rd_TADetectDiscrimSSVEF2(exptDir, sessionDir, fileBase, ...
                        analStr, ssvefFreq, [], iqrThresh, weightChannels, trialSelection);
                    close all;
                end
            else
                error('set either nTopChannels or iqrThresh to a value for channel selection')
            end
        end
        
    end
end
fprintf('done.\n')
    
    
    
    
    
    
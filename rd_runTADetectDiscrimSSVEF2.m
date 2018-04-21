% rd_runTADetectDiscrimSSVEF2.m

% works with SSVEF2, 3, 5, 6, dva (just change number in function call)

%% setup
exptDir = '/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG';
analStr = 'ebi'; % '', 'ebi', etc.
ssvefFreqs = 30; %[30 40];
exptType = 'TADetectDiscrim';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% channel selection: choose nTopChs, iqrThreshs, or wholeBrain
nTopChs = 5; % 1, 5, [1 5], etc.
iqrThreshs = [];
wholeBrain = 0;
weightChannels = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'all','correct','incorrect','validCorrect','detectHit','detectMiss','detectFA','detectCR','discrimCorrect','discrimIncorrect'
trialSelections = {'all'}; 
respTargetSelection = ''; %'T2Resp';

subjects = {'R0817_20150504', 'R0973_20150727', 'R0974_20150728', ...
    'R0861_20150813', 'R0504_20150805', 'R0983_20150813', ...
    'R0898_20150828', 'R0436_20150904', 'R1018_20151118', ...
    'R1019_20151118','R1021_20151120','R1026_20151211', ...
    'R0852_20151211','R1027_20151216','R1028_20151216',...
    'R1029_20151222'}; % N=16
% subjects = subjects(8:end);
nSubjects = numel(subjects);

%% run analysis
% trial selection
for iTS = 1:numel(trialSelections)
    trialSelection = trialSelections{iTS};
    
    % subject
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
            if wholeBrain
                rd_TADetectDiscrimDVA(exptDir, sessionDir, fileBase, ...
                    analStr, ssvefFreq, trialSelection, respTargetSelection);
                close all;
            else
                if ~isempty(nTopChs) && ~isempty(iqrThreshs)
                    error('set either nTopChs or iqrThreshs to empty')
                else
                    if ~isempty(nTopChs)
                        for nTopChannels = nTopChs
                            rd_TADetectDiscrimSSVEF3(exptDir, sessionDir, fileBase, ...
                                analStr, ssvefFreq, nTopChannels, [], weightChannels, ...
                                trialSelection, respTargetSelection, exptType);
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
    end
end
fprintf('done.\n')
    
    

    
    
    
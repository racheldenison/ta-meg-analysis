% rd_runPlotTADetectDiscrimGroup.m

% measures = {'ts', 'w', 'h', 'tf', 'stf', 'w-single', 'ts-single', 'stf-single'};
% measures = {'ts-single'};
% measures = {'ts', 'w', 'tf', 'stf'};
measures = {'w-single', 'itpc-single'};
nMeasures = numel(measures);

trialSelections = {'_allTrials'}; % '', '_allTrials'
% trialSelections = {'detectHit','detectMiss','detectFA','detectCR'}; 
% trialSelections = {'detectHit','detectMiss','detectFA','detectCR','discrimCorrect','discrimIncorrect','validCorrect'}; 
nTrialSelections = numel(trialSelections);

normalizeOptions = {'none','stim'};
nNormOptions = numel(normalizeOptions);

for iM = 1:nMeasures
    measure = measures{iM};
    
    for iTS = 1:nTrialSelections
        trialSelection = trialSelections{iTS};
        selectionStr = sprintf('topChannels5%s', trialSelection);
        
        for iN = 1:nNormOptions
            normalizeOption = normalizeOptions{iN};
            rd_plotTADetectDiscrimGroup(measure, selectionStr, normalizeOption);
            close all
        end
    end
end
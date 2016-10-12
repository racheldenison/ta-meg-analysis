% rd_runPlotTADetectDiscrimGroup.m

% measures = {'ts', 'w', 'h', 'tf', 'stf', 'w-single', 'ts-single', 'stf-single'};
measures = {'ts-single'};
nMeasures = numel(measures);

trialSelections = {'detectHit','detectMiss','detectFA','detectCR'}; 
% trialSelections = {'detectHit','detectMiss','detectFA','detectCR','discrimCorrect','discrimIncorrect','validCorrect'}; 
nTrialSelections = numel(trialSelections);

for iM = 1:nMeasures
    measure = measures{iM};
    
    for iTS = 1:nTrialSelections
        trialSelection = trialSelections{iTS};
        selectionStr = sprintf('topChannels5_%sTrialsT1Resp', trialSelection);
        
        rd_plotTADetectDiscrimGroup(measure, selectionStr);
        close all
    end
end
% rd_runPlotTADetectDiscrimGroup.m

measures = {'ts', 'w', 'h', 'tf', 'stf'};
% measures = {'w','h'};
nMeasures = numel(measures);

for iM = 1:nMeasures
    measure = measures{iM};
    rd_plotTADetectDiscrimGroup(measure);
end
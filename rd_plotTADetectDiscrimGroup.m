% rd_plotTADetectDiscrimGroup.m

%% Setup
exptDir = '/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG';
analStr = 'ebi_ft'; % '', 'ebi', etc.
ssvefFreq = 30;
nTopChannels = 5; % 1, 5, etc.
measure = 'stf'; % w h tf stf

subjects = {'R0817_20150504', 'R0973_20150727', 'R0974_20150728', ...
    'R0861_20150813', ...
    'R0504_20150805', 'R0983_20150813', 'R0898_20150828','R0436_20150904'};
% subjects = {'R0817_20150504', 'R0973_20150727', ...
%     'R0861_20150813', 'R0898_20150828'};
% subjects = {'R0973_20150727', 'R0974_20150728', ...
%     'R0861_20150813', ...
%     'R0504_20150805', 'R0983_20150813', 'R0898_20150828','R0436_20150904'};

nSubjects = numel(subjects);

tstart = -500; % ms
tstop = 3600; % ms
t = tstart:tstop;
evTimes = [0 500 1500 2100 3100];
eventTimes = evTimes;

%% Get data
for iSubject = 1:nSubjects
    subject = subjects{iSubject};
    
    sessionDir = subject;
    dataDir = sprintf('%s/%s', exptDir, sessionDir);
    matDir = sprintf('%s/mat', dataDir);
    
    analysisFile = dir(sprintf('%s/analysis_*_%s_topChannels%d_%dHz.mat', matDir, analStr, nTopChannels, ssvefFreq));
    
    if numel(analysisFile)==1
        load(sprintf('%s/%s', matDir, analysisFile.name))
    else
        error('too many or too few matching analysis files')
    end
    
    switch measure
        case 'w'
            groupData.amps(:,:,iSubject) = A.wAmps;
            groupData.ampsAtt(:,:,iSubject) = A.wAmpsAtt;
            groupData.ampsPA(:,:,iSubject) = A.wAmpsPA;
        case 'h'
            groupData.amps(:,:,iSubject) = A.hAmps;
            groupData.ampsAtt(:,:,iSubject) = A.hAmpsAtt;
            groupData.ampsPA(:,:,iSubject) = A.hAmpsPA;
        case 'tf'
            groupData.amps(:,:,:,iSubject) = A.tfAmps;
            groupData.ampsAtt(:,:,:,iSubject) = A.tfAmpsAtt;
            groupData.ampsPA(:,:,:,iSubject) = A.tfAmpsPA;
            groupData.paDiff(:,:,:,iSubject) = A.tfPADiff;
        case 'stf'
            groupData.amps(:,:,:,iSubject) = A.stfAmps;
            groupData.ampsAtt(:,:,:,iSubject) = A.stfAmpsAtt;
            groupData.ampsPA(:,:,:,iSubject) = A.stfAmpsPA;
            groupData.paDiff(:,:,:,iSubject) = A.stfPADiff;
        otherwise
            error('measure not recognized')
    end
end

%% Calculate group mean and ste
fieldNames = fieldnames(groupData);
nFields = numel(fieldNames);
for iF = 1:nFields
    fieldName = fieldNames{iF};
    vals = groupData.(fieldName);
    sdim = numel(size(vals)); % subject dimension
    groupMean.(fieldName) = mean(vals, sdim);
    groupSte.(fieldName) = std(vals, 0, sdim)./sqrt(nSubjects);
end

%% Plot figs
switch measure
    case {'w','h'}
        rd_plotTADetectDiscrimGroupTS(A, measure, ...
            groupData, groupMean, groupSte)
    case {'tf','stf'}
        rd_plotTADetectDiscrimGroupTimeFreq(A, measure, ...
            groupMean.amps, groupMean.ampsAtt, groupMean.ampsPA, groupMean.paDiff)
    otherwise
        error('measure not recognized')
end


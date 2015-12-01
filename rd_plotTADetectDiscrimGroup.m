function rd_plotTADetectDiscrimGroup(measure)

%% Args
if ~exist('measure','var')
    measure = 'w'; % ts w h tf stf
end

%% Setup
exptDir = '/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG';
analStr = 'ebi_ft'; % '', 'ebi', etc.
ssvefFreq = 30;
nTopChannels = 5; % 1, 5, etc.

% subjects = {'R0817_20150504', 'R0973_20150727', 'R0974_20150728', ...
%     'R0861_20150813', 'R0504_20150805', 'R0983_20150813', ...
%     'R0898_20150828', 'R0436_20150904', 'R0988_20150904', ...
%       'R1018_20151118', 'R1019_20151118', 'R1021_20151120'};
subjects = {'R0817_20150504', 'R0973_20150727', 'R0974_20150728', ...
    'R0861_20150813', 'R0504_20150805', 'R0983_20150813', ...
    'R0898_20150828', 'R0436_20150904', 'R1018_20151118', ...
    'R1019_20151118','R1021_20151120'};
% subjects = {'R0817_20150504', 'R0973_20150727', ...
%     'R0861_20150813', 'R0504_20150805', 'R0898_20150828'}; % endo
% subjects = {'R0973_20150727', 'R0861_20150813', 'R0504_20150805', ...
%     'R0983_20150813', 'R0436_20150904'}; % exo

nSubjects = numel(subjects);

saveFigs = 0;
figDir = sprintf('%s/Group/figures/%s', exptDir, analStr);
figStr = sprintf('gN%d_%dHz_topChannels%d', nSubjects, ssvefFreq, nTopChannels);

tstart = -500; % ms
tstop = 3600; % ms
t = tstart:tstop;
eventTimes = [0 500 1500 2100 3100];

normalizeOption = 'amp'; % 'none','amp','commonBaseline'

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
        case 'ts'
            groupData.tsAmps(:,:,:,iSubject) = A.trigMean;
            groupData.fAmps(:,:,:,iSubject) = A.amps;
            groupData.targetPA(:,:,:,iSubject) = A.targetPA;
            groupData.targetPADiff(:,:,iSubject) = A.targetPADiff;
            groupData.targetPADiffAmps(:,:,iSubject) = A.targetPADiffAmps;
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

%% Normalize to baseline for each subject - common baseline across conds
switch normalizeOption
    case 'commonBaseline'
        bwin = [-500 0];
        fieldNames = fieldnames(groupData);
        nFields = numel(fieldNames);
        for iF = 1:nFields
            fieldName = fieldNames{iF};
            vals = groupData.(fieldName);
            sz = size(vals);
            tdim = find(sz==numel(t));
            if tdim==1
                baseline = nanmean(nanmean(vals(t>=bwin(1) & t<bwin(2),:,:),1),2);
            else
                baseline = nanmean(nanmean(vals(:,t>=bwin(1) & t<bwin(2),:),1),2);
            end
            baselineVals = repmat(baseline, sz(1), sz(2));
            groupData.(fieldName) = (vals - baselineVals)./baselineVals; % change
        end
    case 'amp'
        bwin = [500 3100];
        vals = groupData.amps;
        baselineStim = nanmean(nanmean(vals(t>=bwin(1) & t<bwin(2),1:end-1,:),2),1);
        baselineBlank = nanmean(vals(t>=bwin(1) & t<bwin(2),end,:),1);
        baseline = baselineStim-baselineBlank;
        fieldNames = fieldnames(groupData);
        nFields = numel(fieldNames);
        for iF = 1:nFields
            fieldName = fieldNames{iF};
            vals = groupData.(fieldName);
            sz = size(vals);
            blankVals = repmat(baselineBlank, sz(1), sz(2));
            baselineVals = repmat(baseline, sz(1), sz(2));
            groupData.(fieldName) = (vals - blankVals)./baselineVals; % relative to average amplitude
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

%% A bit of stats
% [h p] = ttest(squeeze(groupData.ampsAtt(1,:,:))', squeeze(groupData.ampsAtt(2,:,:))');
% hold on
% plot(t,h)
% plot(t,-log10(p))
% plot(t, ones(size(t))*(-log10(.05)))

%% Plot figs
switch measure
    case 'ts'
        rd_plotTADetectDiscrimGroupTS(A, ...
            groupMean, ...
            saveFigs, figDir, figStr)
    case {'w','h'}
        rd_plotTADetectDiscrimGroupAmps(A, measure, subjects, ...
            groupData, groupMean, groupSte, ...
            saveFigs, figDir, figStr)
    case {'tf','stf'}
        rd_plotTADetectDiscrimGroupTimeFreq(A, measure, ...
            groupMean, ...
            saveFigs, figDir, figStr)
    otherwise
        error('measure not recognized')
end


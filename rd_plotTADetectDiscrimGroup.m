function rd_plotTADetectDiscrimGroup(measure)

% Args
if ~exist('measure','var')
    measure = 'ts-single'; % ts w h tf stf w-single stf-single
end

% Setup
exptDir = '/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG';
analStr = 'ebi_ft'; % '', 'ebi', etc.
ssvefFreq = 30;
selectionStr = 'topChannels5_allTrials'; %'topChannels5_allTrials'; %'topChannels5'; %'topChannels5_detectHitTrials'; %'topChannels10W_allTrials'; %'topChannels5_validCorrectTrials'; %'iqrThresh10_allTrials';
if strfind(measure,'single')
    aggStr = '_singleTrials';
else
    aggStr = '';
end

subjects = {'R0817_20150504', 'R0973_20150727', 'R0974_20150728', ...
    'R0861_20150813', 'R0504_20150805', 'R0983_20150813', ...
    'R0898_20150828', 'R0436_20150904', 'R1018_20151118', ...
    'R1019_20151118','R1021_20151120','R1026_20151211', ...
    'R0852_20151211','R1027_20151216','R1028_20151216',...
    'R1029_20151222'}; % N=16

% subjects = {'R0817_20150504', 'R0973_20150727', 'R0974_20150728', ...
%     'R0861_20150813', 'R0504_20150805', 'R0983_20150813', ...
%     'R0898_20150828', 'R0436_20150904'}; % N=8 first half

% subjects = {'R1018_20151118', ...
%     'R1019_20151118','R1021_20151120','R1026_20151211', ...
%     'R0852_20151211','R1027_20151216','R1028_20151216',...
%     'R1029_20151222'}; % N=8 second half

% enhancers/suppressers topChannels5_allTrials [-200 200] T1&T2 absent
% subjects = {'R0817_20150504', 'R0973_20150727', ...
%     'R0861_20150813', 'R0504_20150805', ...
%     'R0898_20150828', 'R1018_20151118', ...
%     'R1021_20151120','R1029_20151222'}; % N=8 enhancers (based on absent +/-200 ms)
%
% subjects = {'R0974_20150728', ...
%     'R0983_20150813', ...
%     'R0436_20150904', ...
%     'R1019_20151118','R1026_20151211', ...
%     'R0852_20151211','R1027_20151216','R1028_20151216'}; % N=8 suppressers (based on absent +/-200 ms)

% enhancers/suppressers topChannels10W_allTrials [-200 200] T1&T2 absent
% subjects = {'R0817_20150504', 'R0973_20150727', ...
%     'R0861_20150813', 'R0504_20150805', ...
%     'R0898_20150828', 'R0436_20150904', 'R1018_20151118', ...
%     'R1019_20151118','R1021_20151120','R1026_20151211', ...
%     'R1029_20151222'}; % N=11 enhancers
%
% subjects = {'R0974_20150728', ...
%     'R0983_20150813', ...
%     'R0852_20151211','R1027_20151216','R1028_20151216'}; % N=5 suppressers

% enhancers/suppressers topChannels10W_allTrials [-200 0] T1 absent
% subjects = {'R0817_20150504', 'R0973_20150727', 'R0974_20150728', ...
%     'R0861_20150813', 'R0504_20150805',  ...
%     'R0898_20150828', 'R1018_20151118', ...
%     'R1019_20151118', 'R1029_20151222'}; % N=9 enhancers
%
% subjects = {'R0983_20150813', 'R0436_20150904', ...
%     'R1021_20151120','R1026_20151211', ...
%     'R0852_20151211','R1027_20151216','R1028_20151216'}; % N=7 suppressers


% amp groups (based on topchannels10W all trials)
% subjects = {'R0974_20150728','R0861_20150813', ...
%     'R0983_20150813','R1027_20151216'}; % N=4 low (<100)
% subjects = {'R0817_20150504','R0436_20150904', ...
%     'R1018_20151118','R1019_20151118', 'R1026_20151211', ...
%     'R0852_20151211','R1029_20151222'}; % N=7 med (>100 & <150)
% subjects = {'R0973_20150727','R0504_20150805', ...
%     'R0898_20150828','R1021_20151120','R1028_20151216'}; % N=5 high (>150)


% subjects = {'R0973_20150727','R0983_20150813', ...
%     'R1021_20151120', ...
%     'R0852_20151211','R1027_20151216','R1028_20151216'}; % good behav
% subjects = {'R0817_20150504', 'R0973_20150727', ...
%     'R0861_20150813', 'R0504_20150805', 'R0898_20150828'}; % endo
% subjects = {'R0973_20150727', 'R0861_20150813', 'R0504_20150805', ...
%     'R0983_20150813', 'R0436_20150904'}; % exo

% subjects = {'R0861_20150813','R0504_20150805','R0983_20150813',...
%     'R1021_20151120','R1026_20151211','R0852_20151211',...
%     'R1027_20151216','R1028_20151216','R1029_20151222'}; % discrim1 T1

% subjects = {'R0817_20150504','R0504_20150805','R0983_20150813',...
%     'R0436_20150904','R1018_20151118','R1019_20151118',...
%     'R1021_20151120','R1026_20151211','R1028_20151216',...
%     'R1029_20151222'}; % N=10, subjects with dip (as measured by fitting)

nSubjects = numel(subjects);

saveFigs = 0;
figDir = sprintf('%s/Group/figures/%s', exptDir, analStr);
figStr = sprintf('gN%d%s_%dHz_%s', nSubjects, aggStr, ssvefFreq, selectionStr);

tstart = -500; % ms
tstop = 3600; % ms
t = tstart:tstop;
eventTimes = [0 500 1500 2100 3100];

normalizeOption = 'none'; % 'none','commonBaseline','amp','stim'

%% Get data
for iSubject = 1:nSubjects
    subject = subjects{iSubject};
    
    sessionDir = subject;
    dataDir = sprintf('%s/%s', exptDir, sessionDir);
    matDir = sprintf('%s/mat', dataDir);
    
    analysisFile = dir(sprintf('%s/analysis%s_*_%s_%s_%dHz.mat', matDir, aggStr, analStr, selectionStr, ssvefFreq));

    if isempty(aggStr)
        removeFile = [];
        for i = 1:numel(analysisFile)
            if strfind(analysisFile.name,'singleTrials')
                removeFile = i;
            end
        end
        analysisFile(removeFile) = [];
    end
        
    if numel(analysisFile)==1
        load(sprintf('%s/%s', matDir, analysisFile.name))
    else
        error('%s: too many or too few matching analysis files', subject)
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
        case 'w-single'
            groupData.amps(:,:,iSubject) = squeeze(nanmean(A.wAmps,2));
            groupData.ampsAtt(:,:,iSubject) = squeeze(nanmean(A.wAmpsAtt,2)); % transponse if plotting with amps function
            groupData.ampsPA(:,:,iSubject) = squeeze(nanmean(A.wAmpsPA,2)); % transpose if plotting with amps function
            % comment out if plotting with amps function
            groupData.ampsAll(:,iSubject) = squeeze(nanmean(A.wAmpsAll,2));
            groupData.PAAUT(:,:,:,iSubject) = squeeze(nanmean(A.wPAAUT,2));
            groupData.PAT(:,:,:,iSubject) = squeeze(nanmean(A.wPAT,2));
            groupData.AUT(:,:,:,iSubject) = squeeze(nanmean(A.wAUT,2));
            groupData.PAAU(:,:,iSubject) = squeeze(nanmean(A.wPAAU,2));
            groupData.PA(:,:,iSubject) = squeeze(nanmean(A.wPA,2));
            groupData.AU(:,:,iSubject) = squeeze(nanmean(A.wAU,2));
        case 'stf-single'
            groupData.amps(:,:,:,iSubject) = A.stfAmps;
            groupData.ampsAtt(:,:,:,iSubject) = A.stfAmpsAtt;
            groupData.ampsPA(:,:,:,iSubject) = A.stfAmpsPA;
            groupData.paDiff(:,:,:,iSubject) = A.stfPADiff;
            groupData.PAAUT(:,:,:,:,iSubject) = A.stfPAAUT;
            groupData.PAT(:,:,:,:,iSubject) = A.stfPAT;
            groupData.AUT(:,:,:,:,iSubject) = A.stfAUT;
            groupData.PAAU(:,:,:,iSubject) = A.stfPAAU;
            groupData.PA(:,:,:,iSubject) = A.stfPA;
            groupData.AU(:,:,:,iSubject) = A.stfAU;
        case 'ts-single'
            groupData.tsAmpsS(:,:,:,iSubject) = A.trigMeanMean;
            groupData.paauTS(:,:,:,:,iSubject) = A.paauT;
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
            if numel(sz)==2
                blankVals = repmat(squeeze(baselineBlank)', sz(1), 1);
                baselineVals = repmat(squeeze(baseline)', sz(1), 1);
            elseif numel(sz)==4
                baselineBlank4(1,1,1,1:sz(4)) = squeeze(baselineBlank);
                baseline4(1,1,1,1:sz(4)) = squeeze(baseline);
                blankVals = repmat(baselineBlank4, sz(1), sz(2), sz(3));
                baselineVals = repmat(baseline4, sz(1), sz(2), sz(3));
            else
                blankVals = repmat(baselineBlank, sz(1), sz(2));
                baselineVals = repmat(baseline, sz(1), sz(2));
            end
            groupData.(fieldName) = (vals - blankVals)./baselineVals; % relative to average amplitude
        end
    case 'stim'
        bwin = [500 3100];
        vals = groupData.amps;
        baselineStim = nanmean(nanmean(vals(t>=bwin(1) & t<bwin(2),1:end-1,:),2),1);
        baseline = baselineStim;
        fieldNames = fieldnames(groupData);
        nFields = numel(fieldNames);
        for iF = 1:nFields
            fieldName = fieldNames{iF};
            vals = groupData.(fieldName);
            sz = size(vals);
            baselineVals = repmat(baseline, sz(1), sz(2));
            groupData.(fieldName) = vals./baselineVals; % relative to average stim 
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
%     case {'w','h','w-single'}
    case {'w','h'}
        rd_plotTADetectDiscrimGroupAmps(A, measure, subjects, ...
            groupData, groupMean, groupSte, ...
            saveFigs, figDir, figStr)
    case {'tf','stf','stf-single'}
        rd_plotTADetectDiscrimGroupTimeFreq(A, measure, ...
            groupMean, saveFigs, figDir, figStr)
        if strcmp(measure, 'stf-single') % extra plots
            rd_plotTADetectDiscrimGroupTimeFreqSingle(A, measure, ...
                groupMean, saveFigs, figDir, figStr)
        end
    case 'w-single'
        rd_plotTADetectDiscrimGroupAmpsSingle(A, measure, subjects, ...
            groupData, groupMean, groupSte, ...
            saveFigs, figDir, figStr)
    case 'ts-single'
        rd_plotTADetectDiscrimGroupTSSingle(A, measure, subjects, ...
            groupData, groupMean, groupSte, ...
            saveFigs, figDir, figStr)
    otherwise
        error('measure not recognized')
end


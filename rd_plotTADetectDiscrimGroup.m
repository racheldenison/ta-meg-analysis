function [groupData, groupMean, groupSte, A] = rd_plotTADetectDiscrimGroup(measure, selectionStr, normalizeOption)

%% Args
if ~exist('measure','var') || isempty(measure)
    measure = 'ts-single'; % ts w h tf stf w-single stf-single ts-single itpc-single w-single-wb stf-single-wb itpc-single-wb pre-single-wb
end
if ~exist('selectionStr','var') || isempty(selectionStr)
    selectionStr = 'topChannels5_allTrials'; %'topChannels5'; %'topChannels5_allTrials'; %'wholebrain_allTrials' %'topChannels5_detectHitTrialsT1Resp'; %'topChannels5_allTrials'; %'topChannels5'; %'topChannels5_detectHitTrials'; %'topChannels10W_allTrials'; %'topChannels5_validCorrectTrials'; %'iqrThresh10_allTrials';
end
if ~exist('normalizeOption','var') || isempty(normalizeOption)
    normalizeOption = 'none'; % 'none','commonBaseline','amp','stim'
end

%% Setup
exptType = 'TADetectDiscrim'; % 'TADetectDiscrim','TANoise';

switch exptType
    case 'TADetectDiscrim'
        exptDir = '/Local/Users/denison/Data/TADetectDiscrim/MEG';
%         exptDir = '/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG';
        ssvefFreq = 30;
        t = -500:3600;
        tidx = 1:numel(t);
    case 'TANoise'
        exptDir = '/Local/Users/denison/Data/TANoise/MEG';
        ssvefFreq = 20;
        t = -1500:5700;
        tidx = 1:6701;
        tftidx = 1:671;
    otherwise
        error('exptType not recognized')
end
analStr = 'ebi_ft'; % '', 'ebi', etc.
ssvefStr = sprintf('%dHz',ssvefFreq);
% ssvefStr = 'FOI';

plotFigs = 1;
saveFigs = 0;

switch exptType
    case 'TADetectDiscrim'
        subjects = {'R0817_20150504', 'R0973_20150727', 'R0974_20150728', ...
            'R0861_20150813', 'R0504_20150805', 'R0983_20150813', ...
            'R0898_20150828', 'R0436_20150904', 'R1018_20151118', ...
            'R1019_20151118','R1021_20151120','R1026_20151211', ...
            'R0852_20151211','R1027_20151216','R1028_20151216',...
            'R1029_20151222'}; % N=16 TADetectDiscrim
    case 'TANoise'
        subjects = {'R0817_20171212','R0817_20171213',...
            'R1187_20180105','R1187_20180108',...
            'R0983_20180111','R0983_20180112',...
            'R0898_20180112','R0898_20180116',...
            'R1021_20180208','R1021_20180212',...
            'R1103_20180213','R1103_20180215',...
            'R0959_20180219','R0959_20180306'}; % N=7 x 2 sessions TANoise
    otherwise
        error('exptType not recognized')
end

% subjects = subjects([1:4 7:14]);
% subjects = subjects([1 2 4 5 7 8 10 12 14 16]);

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

% enhancers/suppressers topChannels5_allTrials_normStim [-250 -50] T1&T2 absent
% subjects = {'R0817_20150504', 'R0973_20150727', ...
%     'R0861_20150813', 'R0504_20150805', ...
%     'R0898_20150828', 'R1018_20151118', ...
%     'R1021_20151120','R0852_20151211','R1029_20151222'}; % N=9 enhancers 
% 
% subjects = {'R0974_20150728', ...
%     'R0983_20150813', ...
%     'R0436_20150904', ...
%     'R1019_20151118','R1026_20151211', ...
%     'R1027_20151216','R1028_20151216'}; % N=7 suppressers 

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

% subjects = {'R0817_20150504', 'R0973_20150727', ...
%     'R0504_20150805', ...
%     'R0898_20150828', 'R0436_20150904', 'R1018_20151118', ...
%     'R1019_20151118','R1021_20151120','R1026_20151211', ...
%     'R0852_20151211','R1028_20151216',...
%     'R1029_20151222'}; % N=12, subjects with reasonable signal

nSubjects = numel(subjects);

if strfind(measure,'single')
    aggStr = '_singleTrials';
else
    aggStr = '';
end
normStr = sprintf('_norm%s',[upper(normalizeOption(1)) normalizeOption(2:end)]);

figDir = sprintf('%s/Group/figures/%s', exptDir, analStr);
figStr = sprintf('gN%d%s_%s_%s%s', nSubjects, aggStr, ssvefStr, selectionStr, normStr);

%% Get data
for iSubject = 1:nSubjects
    subject = subjects{iSubject};
    
    sessionDir = subject;
    dataDir = sprintf('%s/%s', exptDir, sessionDir);
    matDir = sprintf('%s/mat', dataDir);
    
    analysisFile = dir(sprintf('%s/analysis%s_*_%s_%s_%s.mat', matDir, aggStr, analStr, selectionStr, ssvefStr));

    if isempty(aggStr)
        removeFile = [];
        for i = 1:numel(analysisFile)
            if strfind(analysisFile(i).name,'singleTrials')
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
            groupData.amps(:,:,iSubject) = A.wAmps(tidx,:);
            groupData.ampsAtt(:,:,iSubject) = A.wAmpsAtt(:,tidx);
            groupData.ampsPA(:,:,iSubject) = A.wAmpsPA(:,tidx);
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
            groupData.amps(:,:,:,iSubject) = A.stfAmps(:,tftidx,:);
            groupData.ampsAtt(:,:,:,iSubject) = A.stfAmpsAtt(:,tftidx,:);
            groupData.ampsPA(:,:,:,iSubject) = A.stfAmpsPA(:,tftidx,:);
            groupData.paDiff(:,:,:,iSubject) = A.stfPADiff(:,tftidx,:);
        case 'w-single'
            groupData.amps(:,:,iSubject) = squeeze(nanmean(A.wAmps(tidx,:,:),2));
            groupData.ampsAtt(:,:,iSubject) = squeeze(nanmean(A.wAmpsAtt(tidx,:,:),2))'; % transponse if plotting with amps function
            groupData.ampsPA(:,:,iSubject) = squeeze(nanmean(A.wAmpsPA(tidx,:,:),2))'; % transpose if plotting with amps function
            % comment out if plotting with amps function
%             groupData.ampsAll(:,iSubject) = squeeze(nanmean(A.wAmpsAll,2));
%             groupData.PAAUT(:,:,:,iSubject) = squeeze(nanmean(A.wPAAUT,2));
%             groupData.PAT(:,:,:,iSubject) = squeeze(nanmean(A.wPAT,2));
%             groupData.AUT(:,:,:,iSubject) = squeeze(nanmean(A.wAUT,2));
%             groupData.PAAU(:,:,iSubject) = squeeze(nanmean(A.wPAAU,2));
%             groupData.PA(:,:,iSubject) = squeeze(nanmean(A.wPA,2));
%             groupData.AU(:,:,iSubject) = squeeze(nanmean(A.wAU,2));
        case 'itpc-single'
            groupData.amps(:,:,iSubject) = A.wITPC(tidx,:);
            groupData.ampsAtt(:,:,iSubject) = A.wITPCAtt(tidx,:)';
            groupData.ampsPA(:,:,iSubject) = A.wITPCPA(tidx,:)';
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
            groupData.tsAmps(:,:,:,iSubject) = squeeze(nanmean(A.trigMean,3));
            groupData.fAmps(:,:,:,iSubject) = squeeze(nanmean(A.amps,3));
            groupData.tsAmpsS(:,:,:,iSubject) = A.trigMeanMean;
            groupData.paauTS(:,:,:,:,iSubject) = A.paauT;
            % number of trials
            excludedTrials = squeeze(all(isnan(A.trigMeanMean),1));
            groupData.nTrialsPerCond(:,iSubject) = sum(1-excludedTrials);
        case 'w-single-wb'
            groupData.amps(:,:,:,iSubject) = A.wAmps;
            groupData.ampsAll(:,:,iSubject) = A.wAmpsAll;
            groupData.PAAUT(:,:,:,:,iSubject) = A.wPAAUT;
            groupData.PAT(:,:,:,:,iSubject) = A.wPAT;
            groupData.AUT(:,:,:,:,iSubject) = A.wAUT;
            groupData.PAAU(:,:,:,iSubject) = A.wPAAU;
            groupData.PA(:,:,:,iSubject) = A.wPA;
            groupData.AU(:,:,:,iSubject) = A.wAU;
        case 'stf-single-wb'
            groupData.amps(:,:,:,:,iSubject) = A.stfAmps;
            groupData.ampsAtt(:,:,:,:,iSubject) = A.stfAmpsAtt;
            groupData.ampsPA(:,:,:,:,iSubject) = A.stfAmpsPA;
            groupData.paDiff(:,:,:,:,iSubject) = A.stfPADiff;
            groupData.PAAUT(:,:,:,:,:,iSubject) = A.stfPAAUT;
            groupData.PAT(:,:,:,:,:,iSubject) = A.stfPAT;
            groupData.AUT(:,:,:,:,:,iSubject) = A.stfAUT;
            groupData.PAAU(:,:,:,:,iSubject) = A.stfPAAU;
            groupData.PA(:,:,:,:,iSubject) = A.stfPA;
            groupData.AU(:,:,:,:,iSubject) = A.stfAU;  
        case 'itpc-single-wb'
            groupData.amps(:,:,:,:,iSubject) = A.itpcAmps;
            groupData.PAAUT(:,:,:,:,:,iSubject) = A.itpcPAAUT;
            groupData.PAT(:,:,:,:,:,iSubject) = A.itpcPAT;
            groupData.AUT(:,:,:,:,:,iSubject) = A.itpcAUT;
            groupData.PAAU(:,:,:,:,iSubject) = A.itpcPAAU;
            groupData.PA(:,:,:,:,iSubject) = A.itpcPA;
            groupData.AU(:,:,:,:,iSubject) = A.itpcAU;   
        case 'pre-single-wb'
            groupData.alpha(:,:,:,iSubject) = A.alpha.preAUTZ;
            groupData.ssvef30(:,:,:,iSubject) = A.ssvef30.preAUTZ;
            groupData.ssvef40(:,:,:,iSubject) = A.ssvef40.preAUTZ;
            groupData.broadband(:,:,:,iSubject) = A.broadband.preAUTZ;
        otherwise
            error('measure not recognized')
    end
end

switch measure
    case {'stf-single-wb','itpc-single-wb'}
        groupDataDiff.AU = squeeze(groupData.AU(:,:,:,1,:)-groupData.AU(:,:,:,2,:));
        groupDataDiff.PA = squeeze(groupData.PA(:,:,:,1,:)-groupData.PA(:,:,:,2,:));
        groupDataDiff.AUT = squeeze(groupData.AUT(:,:,:,1,:,:)-groupData.AUT(:,:,:,2,:,:));
        dotstat = 1;
    case 'pre-single-wb'
        fieldNames = fields(groupData);
        for iF = 1:numel(fieldNames)
            groupDataDiff.(fieldNames{iF}) = squeeze(groupData.(fieldNames{iF})(:,1,:,:) - ...
                groupData.(fieldNames{iF})(:,2,:,:));
        end
        dotstat = 1;
    otherwise
        dotstat = 0;
end

%% Normalize to baseline for each subject - common baseline across conds
A.normalizeOption = normalizeOption;
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
        if strfind(measure,'wb')
            baselineStim = nanmean(nanmean(vals(t>=bwin(1) & t<bwin(2),:,1:end-1,:),3),1);
        else
            baselineStim = nanmean(nanmean(vals(t>=bwin(1) & t<bwin(2),1:end-1,:),2),1);
        end
        baseline = baselineStim;
        fieldNames = fieldnames(groupData);
        nFields = numel(fieldNames);
        for iF = 1:nFields
            fieldName = fieldNames{iF};
            vals = groupData.(fieldName);
            sz = size(vals);
            if strfind(measure,'wb')
                if numel(sz)==3
                    baselineVals = squeeze(repmat(baseline, sz(1), 1, 1, 1));
                elseif numel(sz)==4
                    baselineVals = repmat(baseline, sz(1), 1, sz(3), 1);
                elseif numel(sz)==5
                    baseline5(1,:,1,1,:) = baseline;
                    baselineVals = repmat(baseline5, sz(1), 1, sz(3), sz(4), 1);
                else
                    error('baseline not currently supported for matrix dim ')
                end
            else
                % baselineVals = repmat(baseline, sz(1), sz(2));
                if numel(sz)==2
                    baselineVals = repmat(squeeze(baseline)', sz(1), 1);
                elseif numel(sz)==4
                    baseline4(1,1,1,1:sz(4)) = squeeze(baseline);
                    baselineVals = repmat(baseline4, sz(1), sz(2), sz(3));
                elseif numel(sz)==5
                    baseline5(1,1,1,1,1:sz(5)) = squeeze(baseline);
                    baselineVals = repmat(baseline5, sz(1), sz(2), sz(3), sz(4));
                elseif numel(sz)>5
                    error('baseline for matrix dim>5 not currently supported')
                else
                    baselineVals = repmat(baseline, sz(1), sz(2));
                end
            end
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

%% Calculate group t-stat (uncorrected)
if dotstat
    fieldNames = fieldnames(groupDataDiff);
    nFields = numel(fieldNames);
    for iF = 1:nFields
        fieldName = fieldNames{iF};
        vals = groupDataDiff.(fieldName);
        sdim = numel(size(vals)); % subject dimension
        [h p ci stat] = ttest(vals,0,'dim',sdim);
        groupTStat.(fieldName) = stat.tstat;
    end
end

%% A bit of stats
% [h p] = ttest(squeeze(groupData.ampsAtt(1,:,:))', squeeze(groupData.ampsAtt(2,:,:))');
% hold on
% plot(t,h)
% plot(t,-log10(p))
% plot(t, ones(size(t))*(-log10(.05)))

%% Adjust t if neeeded
A.t = A.t(tidx);

%% Plot figs
if plotFigs
switch measure
    case 'ts'
        rd_plotTADetectDiscrimGroupTS(A, measure, groupMean, ...
            saveFigs, figDir, figStr)
%     case {'w','h','w-single'}
    case {'w','h','itpc-single','w-single'}
        rd_plotTADetectDiscrimGroupAmps2(A, measure, subjects, ...
            groupData, groupMean, groupSte, ...
            saveFigs, figDir, figStr)
    case {'tf','stf','stf-single'}
        rd_plotTADetectDiscrimGroupTimeFreq(A, measure, ...
            groupMean, saveFigs, figDir, figStr)
        if strcmp(measure, 'stf-single') % extra plots
            rd_plotTADetectDiscrimGroupTimeFreqSingle(A, measure, ...
                groupMean, saveFigs, figDir, figStr)
        end
    case 'w-single----'
        rd_plotTADetectDiscrimGroupAmpsSingle(A, measure, subjects, ...
            groupData, groupMean, groupSte, ...
            saveFigs, figDir, figStr)
    case 'ts-single'
        rd_plotTADetectDiscrimGroupTSSingle(A, measure, subjects, ...
            groupData, groupMean, groupSte, ...
            saveFigs, figDir, figStr, selectionStr)
    case 'w-single-wb'
        rd_plotTADetectDiscrimGroupAmpsWholebrain(A, measure, subjects, ...
            groupData, groupMean, groupSte, ...
            saveFigs, figDir, figStr)
    case {'stf-single-wb','itpc-single-wb'}
        rd_plotTADetectDiscrimGroupTimeFreqWholebrain(A, measure, subjects, ...
            groupData, groupMean, groupSte, groupTStat, ...
            saveFigs, figDir, figStr)
    case 'pre-single-wb'
        rd_plotTADetectDiscrimGroupPreWholebrain(A, measure, subjects, ...
            groupData, groupMean, groupSte, groupTStat, ...
            saveFigs, figDir, figStr)
    otherwise
        error('measure not recognized')
end
end


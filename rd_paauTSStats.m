% rd_paauTSStats.m

%% load data
% load('/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG/Group/mat/paauTS_stats_workspace_20160924.mat')
load('/Local/Users/denison/Data/TADetectDiscrim/MEG/Group/mat/paauTS_stats_workspace_20160924.mat')

%% setup
twin = A.targetWindow;
twindow = twin(1):twin(end);

toi = [0 twin(end)];
twoi = find(twindow==toi(1)):find(twindow==toi(end));
nt = numel(twoi);

nSubjects = numel(subjects);
nConds = numel(condNames);
nShuffles = 5000;

%% generate null distribution for full ANOVA F values
% vals is time x conds x subjects
tic
for iShuffle = 1:nShuffles
    fprintf('shuffle %d \t%s\n', iShuffle, datestr(now))
    % shuffle subject means independently
    for iS = 1:nSubjects
        idx = randperm(nConds);
        valsShuffle(:,:,iS) = vals(twoi,idx,iS);
    end

    for it = 1:nt
        data = squeeze(valsShuffle(it,:,:))'; % subjects x conds
        [fvalsNull(it,:,iShuffle), pvalsNull(it,:,iShuffle)] = rd_rmANOVA(data, condNames, factorNames, nLevels);
    end
end
toc

%% generate null distribution for main effect differences
% empirical
paDiff = squeeze(mean(mean(groupDataB.paDiff,2),3)); % time x subjects
auDiff = squeeze(mean(groupDataB.auDiff,2));

% only include after baseline
paDiffMean = mean(paDiff(twoi,:),2);
auDiffMean = mean(auDiff(twoi,:),2);

%% stage 1: calculate CIs across time
for iShuffle = 1:nShuffles
    % shuffle subject means independently
    for iS = 1:nSubjects
        idx = randperm(nConds);
        valsShuffle(:,:,iS) = vals(twoi,idx,iS);
    end

    paDiffSh = squeeze(mean(valsShuffle(:,[1:2 5:6],:),2) - mean(valsShuffle(:,[3:4 7:8],:),2));
    auDiffSh = squeeze(mean(valsShuffle(:,[1 3 5 7],:),2) - mean(valsShuffle(:,[2 4 6 8],:),2));
    
    paDiffNull(:,iShuffle) = mean(paDiffSh,2);
    auDiffNull(:,iShuffle) = mean(auDiffSh,2);
end

paDiffCI = prctile(paDiffNull,[2.5 97.5],2);
auDiffCI = prctile(auDiffNull,[2.5 97.5],2);

paThresh = paDiffMean>paDiffCI(:,2) | paDiffMean<paDiffCI(:,1);
auThresh = auDiffMean>auDiffCI(:,2) | auDiffMean<auDiffCI(:,1);

% empirical cluster sum
[~, paDiffCluster] = rd_clusterSum(paDiffMean, paThresh);
[~, auDiffCluster] = rd_clusterSum(auDiffMean, auThresh);


figure
hold on
plot(auDiffNull)
plot(auDiffMean,'k','LineWidth',2)
plot(auDiffCI,'LineWidth',2)
plot(5*auThresh,'k','LineWidth',2)

figure
hold on
plot(paDiffNull)
plot(paDiffMean,'k','LineWidth',2)
plot(paDiffCI,'LineWidth',2)
plot(5*paThresh,'k','LineWidth',2)

%% stage 2: max cluster sum
% use CIs calculated in previous stage
for iShuffle = 1:nShuffles
    % shuffle subject means independently
    for iS = 1:nSubjects
        idx = randperm(nConds);
        valsShuffle(:,:,iS) = vals(twoi,idx,iS);
    end

    paDiffShMean = mean(squeeze(mean(valsShuffle(:,[1:2 5:6],:),2) - mean(valsShuffle(:,[3:4 7:8],:),2)),2);
    auDiffShMean = mean(squeeze(mean(valsShuffle(:,[1 3 5 7],:),2) - mean(valsShuffle(:,[2 4 6 8],:),2)),2);
    
    paThreshSh = paDiffShMean>paDiffCI(:,2) | paDiffShMean<paDiffCI(:,1);
    auThreshSh = auDiffShMean>auDiffCI(:,2) | auDiffShMean<auDiffCI(:,1);
    
    [~, paDiffClusterNull(iShuffle)] = rd_clusterSum(paDiffShMean, paThreshSh);
    [~, auDiffClusterNull(iShuffle)] = rd_clusterSum(auDiffShMean, auThreshSh);
end

paDiffClusterCI = prctile(paDiffClusterNull,[2.5 97.5],2);
auDiffClusterCI = prctile(auDiffClusterNull,[2.5 97.5],2);

%% max cluster sum of t-stat
[h p ci stats] = ttest(paDiff(twoi,:)');
paDiffTStat = stats.tstat;

[h p ci stats] = ttest(auDiff(twoi,:)');
auDiffTStat = stats.tstat;

% t threshold
tthresh = abs(tinv(.05/2,nSubjects-1));

% empirical cluster sum
[~, paDiffCluster] = rd_clusterSum(paDiffTStat, abs(paDiffTStat)>tthresh);
[~, auDiffCluster] = rd_clusterSum(auDiffTStat, abs(auDiffTStat)>tthresh);

for iShuffle = 1:nShuffles
    % shuffle subject means independently
    for iS = 1:nSubjects
        idx = randperm(nConds);
        valsShuffle(:,:,iS) = vals(twoi,idx,iS);
    end

    paDiffShData = squeeze(mean(valsShuffle(:,[1:2 5:6],:),2) - mean(valsShuffle(:,[3:4 7:8],:),2));
    auDiffShData = squeeze(mean(valsShuffle(:,[1 3 5 7],:),2) - mean(valsShuffle(:,[2 4 6 8],:),2));
    
    [h p ci stats] = ttest(paDiffShData');
    paDiffShTStat = stats.tstat;
    [h p ci stats] = ttest(auDiffShData');
    auDiffShTStat = stats.tstat;
    audshts(:,iShuffle) = auDiffShTStat;
    
    [~, paDiffClusterNull(iShuffle)] = rd_clusterSum(paDiffShTStat, abs(paDiffShTStat)>tthresh);
    [~, auDiffClusterNull(iShuffle)] = rd_clusterSum(auDiffShTStat, abs(auDiffShTStat)>tthresh);
end

paDiffClusterCI = prctile(paDiffClusterNull,95,2);
auDiffClusterCI = prctile(auDiffClusterNull,95,2);






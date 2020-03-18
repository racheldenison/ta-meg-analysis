% rd_itpcStats.m

%% load data
exptDir = pathToTANoise('MEG');

dataFile = sprintf('%s/Group/mat/gN10_itpcAtt_20Hz.mat', exptDir);

D = load(dataFile);

%% setup
t = D.t;
data = D.itpc;
nConds = size(data,2);
nSubjects = size(data,3);

dataDiff = squeeze(data(:,1,:) - data(:,2,:));

%% max cluster sum of t-stat
[h p ci stats] = ttest(dataDiff');
diffTStat = stats.tstat;

% t threshold
tthresh = abs(tinv(.05/2,nSubjects-1));

% empirical cluster sum
[clusters, diffCluster, C] = rd_clusterSum(diffTStat, abs(diffTStat)>tthresh);

%% shuffles
nShuffles = 1000;

for iShuffle = 1:nShuffles
    % shuffle subject means independently
    dataSh = [];
    for iS = 1:nSubjects
        idx = randperm(nConds);
        dataSh(:,:,iS) = data(:,idx,iS);
    end

    dataDiffSh = squeeze(dataSh(:,1,:) - dataSh(:,2,:));
    
    [h p ci stats] = ttest(dataDiffSh');
    diffShTStat = stats.tstat;
    
    [~, diffClusterNull(iShuffle)] = rd_clusterSum(diffShTStat, abs(diffShTStat)>tthresh);
end

diffClusterCI = prctile(diffClusterNull,95,2);

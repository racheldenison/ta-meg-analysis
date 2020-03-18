% rd_TAANOVAStats.m

%% setup
exptDir = pathToTANoise('MEG');
analysisDir = sprintf('%s/Group/mat', exptDir);
analysisName = 'TANoise_N10_itpc_40Hz';
statName = 'Fval';
threshName = 'pval';
threshVal = 0.05;

%% empirical
statTable = readtable(sprintf('%s/%s_ANOVA_%s.txt', analysisDir, analysisName, statName), 'ReadRowNames',false);
threshTable = readtable(sprintf('%s/%s_ANOVA_%s.txt', analysisDir, analysisName, threshName), 'ReadRowNames',false);

t = statTable.time;

names = statTable.Properties.VariableNames;
nTests = numel(names);

clusterStatsEmp = [];
maxAbsClusterStatEmp = [];
CEmp = [];
for iTest = 2:nTests
    name = names{iTest};
    [clusterStatsEmp.(name), maxAbsClusterStatEmp.(name), CEmp.(name)] = rd_clusterStat(statTable.(name), threshTable.(name)<threshVal);
end

%% null
statTable = readtable(sprintf('%s/%s_permutation_ANOVA_%s.txt', analysisDir, analysisName, statName),...
    'ReadRowNames',false);
threshTable = readtable(sprintf('%s/%s_permutation_ANOVA_%s.txt', analysisDir, analysisName, threshName),...
    'ReadRowNames',false);

names = statTable.Properties.VariableNames;
nTests = numel(names);

nPerm = numel(unique(statTable.permutation));

maxAbsClusterStatNull = [];
for iP = 1:nPerm
    if mod(iP,100)==0
        fprintf('%d\n',iP)
    end
    w = statTable.permutation==iP;
    for iTest = 3:nTests
        name = names{iTest};
        [~, maxAbsClusterStatNull.(name)(iP,1)] = rd_clusterStat(statTable.(name)(w,:), threshTable.(name)(w,:)<threshVal);
    end
end

%% cluster p-value
clusterPVal = [];
for iTest = 3:nTests
    name = names{iTest};
    clusterPVal.(name) = nnz(maxAbsClusterStatNull.(name)>maxAbsClusterStatEmp.(name))/nPerm;
end

%% plot
figure
for iTest = 3:nTests
    name = names{iTest};
    subplot(nTests-2,1,iTest-2)
    histogram(maxAbsClusterStatNull.(name))
    vline(maxAbsClusterStatEmp.(name),'r');
    title(und2space(name))
end



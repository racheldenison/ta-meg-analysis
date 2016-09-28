function [clusterSums, maxAbsClusterSum] = rd_clusterSum(ts, thresh)

% example inputs
% ts = rand(1,100);
% thresh = ts>.8;

clusterStarts = find(diff(thresh)==1)+1;
clusterEnds = find(diff(thresh)==-1);

sz = size(clusterStarts);
if sz(1)>sz(2)
    clusterStarts = clusterStarts';
    clusterEnds = clusterEnds';
end

if thresh(1)==1
    clusterStarts = [1 clusterStarts];
end
if thresh(end)==1
    clusterEnds = [clusterEnds length(ts)];
end

if numel(clusterStarts)~=numel(clusterEnds)
    error('cluster starts and ends not equal in number')
end

nClusters = numel(clusterStarts);

clusterSums = 0;
for iC = 1:nClusters
    vals = ts(clusterStarts(iC):clusterEnds(iC));
    clusterSums(iC) = sum(vals);
end

maxAbsClusterSum = max(abs(clusterSums));
function [clusterStats, maxAbsClusterStat, C] = rd_clusterStat(im, thresh, stat)
%
% function [clusterStats, maxAbsClusterStat, C] = rd_clusterStat(im, thresh, [stat])
%
% example inputs
% ts = rand(1,100);
% thresh = ts>.8;
% stat = 'sum';

if nargin < 3
    stat = 'sum';
end

C = bwconncomp(thresh);

clusterStats = 0;
for iC = 1:C.NumObjects
    vals = im(C.PixelIdxList{iC});
    switch stat
        case 'sum'
            clusterStats(iC) = sum(vals);
        case 'mean'
            clusterStats(iC) = mean(vals);
        otherwise
            error('stat not recognized')
    end
end

maxAbsClusterStat = max(abs(clusterStats));
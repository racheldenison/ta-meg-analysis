function dataNorm = normalizeDC(data)

% function dataNorm = normalizeDC(data)
%
% data must be 2 or 3 dimensional with subjects as last dimension

dim = numel(size(data));

switch dim
    case 2
        % data is conditions x subjects
        nConds = size(data,1);
        subjectMean = nanmean(data,1);
        grandMean = nanmean(nanmean(data));
        dc = repmat(grandMean - subjectMean,nConds,1);
        dataNorm = data + dc;
        
    case 3
        % data is conditions1 x conditions2 x subjects
        nConds1 = size(data,1);
        nConds2 = size(data,2);
        subjectMean = nanmean(nanmean(data,1),2);
        grandMean = nanmean(data(:));
        dc = repmat(grandMean - subjectMean, [nConds1, nConds2, 1]);
        dataNorm = data + dc;
        
    otherwise
        error('data must be 2- or 3-dimensional')
end
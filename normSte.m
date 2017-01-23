function ste = normSte(data)

% function ste = normSte(data)
%
% subjects must be in the last dimension


dim = numel(size(data));
nSubjects = size(data,dim);

dataNorm = normalizeDC(data);

ste = std(dataNorm,0,dim)./sqrt(nSubjects);

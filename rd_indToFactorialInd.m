function varargout = rd_indToFactorialInd(ind, levels)

factmat = fullfact(levels);
factInds = factmat(ind,:);

varargout = num2cell(factInds);
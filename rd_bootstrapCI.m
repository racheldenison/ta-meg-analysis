function [ci, err] = rd_bootstrapCI(y)
%
% function [ci, err] = rd_bootstrapCI(y)
%
% ci is 95% confidence interval
% err is distances to upper and lower error bars

fun = @nanmean;
confbounds = [2.5 97.5];
nboot = 100;

m = bootstrp(nboot, fun, y);

ci = prctile(m,confbounds);

err(1,:) = ci(1,:)-fun(y);
err(2,:) = fun(y)-ci(2,:);
function [ci, emp, err, m] = rd_bootstrapCI(y, fun)
%
% function [ci, emp, err, m] = rd_bootstrapCI(y, [fun])
%
% Inputs:
% y is the data. should be oriented so that the mean / other statistic can
% be taken across the first dimension.
% fun is an optional function handle for the statistic you want to
% calculate. default is nanmean.
%
% Outputs:
% ci is the 95% confidence interval
% emp is the empirical (true) mean / other statistic of the data
% err is distances to upper and lower error bars
% m is all of the bootstrapped means / other statistics

if nargin==1
    fun = @nanmean;
end
confbounds = [2.5 97.5];
nboot = 100;

m = bootstrp(nboot, fun, y);

ci = prctile(m,confbounds);

emp = fun(y);

err(1,:) = ci(1,:)-emp;
err(2,:) = emp-ci(2,:);


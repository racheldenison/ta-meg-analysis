% rd_plotTransitionTrials.m

%% find transition trials -- flicker trials following blanks
% one way
flicker1Trials = repmat((2:5:41),14,1) + repmat(((0:13)*41)',1,8);
flicker1Trials = sort(flicker1Trials(:));

% % or another
% flicker1Trials = zeros(size(wBlank));
% nTrials = size(wBlank,1);
% for i=2:nTrials
%     if wBlank(i)==0 && wBlank(i-1)==1
%         flicker1Trials(i) = 1;
%     end
% end
% flicker1Trials = logical(flicker1Trials);

%% look at mean tseries
channels = 1;
flicker1Data = squeeze(mean(trigData(:,channels,flicker1Trials),2)); % mean across channels
        
flicker1Mean = nanmean(flicker1Data,2);

%% plot
tsFigPos = [0 500 1250 375];

figure
set(gcf,'Position',tsFigPos)
plot(t, flicker1Mean)
vline(0,'color','k','LineStyle',':')
xlim([-500 1000])
xlabel('time (ms)')
title(sprintf('%s Ch %d', und2space(sessionDir), channels))
turnwhite


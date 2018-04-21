function rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes)

if nargin < 5
    eventTimes = [];
end

set(gca,'YDir','normal')
set(gca,'XTick',xtick)
set(gca,'XTickLabel',toi(xtick))
set(gca,'YTick',ytick)
set(gca,'YTickLabel',foi(ytick))
hold on
for iEv = 1:numel(eventTimes)
    vline(find(abs(toi-eventTimes(iEv)/1000)<1e-4),'k');
end
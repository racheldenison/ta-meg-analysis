% rd_plotGroupTimeSeries.m

d1 = load('R0817_20150504/wAmpsAttT1T2.mat');
d2 = load('R0973_20150727/wAmpsAttT1T2.mat');
d3 = load('R0974_20150728/wAmpsAttT1T2.mat');
d4 = load('R0504_20150805/wAmpsAttT1T2.mat');

t = d2.t;
d2 = rmfield(d2,'t');
d4 = rmfield(d4,'t');

d(1) = d1;
d(2) = d2;
d(3) = d3;
d(4) = d4;

nSubjects = numel(d);
selectedSubjects = [1:2 4];

evTimes = [0 500 1500 2100 3100];

for iS = 1:nSubjects
    wAmps(iS,:,1) = d(iS).wAmpsAttT1;
    wAmps(iS,:,2) = d(iS).wAmpsAttT2;
end

for iS = 1:nSubjects
    wAmpsN(iS,:,:) = wAmps(iS,:,:)./max(max(wAmps(iS,:,:)));
    wAmpsDiff(iS,:) = wAmps(iS,:,1) - wAmps(iS,:,2);
end

%% summary stats
wAmpsMean = squeeze(mean(wAmps(selectedSubjects,:,:)));
wAmpsSte = squeeze(std(wAmps(selectedSubjects,:,:))./sqrt(numel(selectedSubjects)));

wAmpsDiffMean = squeeze(mean(wAmpsDiff(selectedSubjects,:)));
wAmpsDiffSte = squeeze(std(wAmpsDiff(selectedSubjects,:))./sqrt(numel(selectedSubjects)));

%% plot figs
colors = {'b','r'};
xlims = [t(1) t(end)];
ylims = [-.5 2.5];
diffLims = [-1 1];

% indiv subjects
figure
set(gcf,'Position',[30 450 1200 250])
for iS = 1:nSubjects
    subplot(1,nSubjects,iS)
    hold on
    for i=1:2
        plot(t, wAmps(iS,:,i), colors{i})
    end
    plot(xlims, [0 0], 'k')
    xlim(xlims)
    ylim(ylims)
    for iEv = 1:numel(evTimes)
        vline(evTimes(iEv),'color','k','LineStyle',':');
    end
end

% figure
% set(gcf,'Position',[30 450 1200 250])
% for iS = 1:nSubjects
%     subplot(1,nSubjects,iS)
%     hold on
%     for i=1:2
%         plot(t, wAmpsN(iS,:,i), colors{i})
%     end
%     plot(xlims, [0 0], 'k')
%     xlim(xlims)
%     ylim(ylims)
%     for iEv = 1:numel(evTimes)
%         vline(evTimes(iEv),'color','k','LineStyle',':');
%     end
% end

figure
set(gcf,'Position',[30 450 1200 250])
for iS = 1:nSubjects
    subplot(1,nSubjects,iS)
    hold on
    plot(t, wAmpsDiff(iS,:), 'k', 'LineWidth', 2)
    plot(xlims, [0 0], 'k')
    xlim(xlims)
    ylim(diffLims)
    for iEv = 1:numel(evTimes)
        vline(evTimes(iEv),'color','k','LineStyle',':');
    end
end

% group
figure
hold on
for i=1:2
    shadedErrorBar(t, wAmpsMean(:,i), wAmpsSte(:,i), colors{i}, 1)
end
plot(xlims, [0 0], 'k')
for iEv = 1:numel(evTimes)
    vline(evTimes(iEv),'color','k','LineStyle',':');
end

figure
hold on
for i=1:2
    shadedErrorBar(t, wAmpsMean(:,i), wAmpsDiffSte, colors{i}, 1)
end
plot(xlims, [0 0], 'k')
for iEv = 1:numel(evTimes)
    vline(evTimes(iEv),'color','k','LineStyle',':');
end

figure
hold on
shadedErrorBar(t, wAmpsDiffMean, wAmpsDiffSte, 'k', 1)
plot(xlims, [0 0], 'k')
ylim([-.5 .5])
for iEv = 1:numel(evTimes)
    vline(evTimes(iEv),'color','k','LineStyle',':');
end


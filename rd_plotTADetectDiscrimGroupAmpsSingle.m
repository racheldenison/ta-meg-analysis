function rd_plotTADetectDiscrimGroupAmpsSingle(A, measure, subjects, groupData, groupMean, groupSte, saveFigs, figDir, figStr)

%% args
if nargin<7
    saveFigs = 0;
end
if nargin<9
    figStr = '';
    if saveFigs==1
        error('If youre saving figs, you should specify a figStr')
    end
end

figTitle = und2space(figStr);

%% setup
plotOrder = [1 5 3 7 2 6 4 8 9];

tsFigPos = [0 500 1250 375];
paau3FigPos = [800 30 450 830];
paau6FigPos = [0 90 1000 800];

eventTimes = A.eventTimes;
trigNames = A.trigNames;
attNames = A.attNames;
paNames = A.PANames;
t = A.t;
twin = A.wtwin;

nTrigs = numel(trigNames);
nSubjects = size(groupData.amps, 3);

% colors
extendedMap = flipud(lbmap(nTrigs-1+4,'RedBlue'));
selectedMap = extendedMap([1:(nTrigs-1)/2 (end-(nTrigs-1)/2)+1:end],:);
trigColors = [selectedMap; 0 0 0];
trigBlue = mean(selectedMap(1:(nTrigs-1)/2,:));
trigRed = mean(selectedMap((end-(nTrigs-1)/2)+1:end,:));
trigColorsAtt2 = [trigBlue; trigRed];
trigColorsPA4 = [.52 .37 .75; .31 .74 .40; .27 .51 .84; 1.0 .57 .22];

set(0,'defaultLineLineWidth',1)

switch measure
    case 'w-single'
        ylims = [300 700];
    otherwise
        error('measure not recognized')
end

%% PAAU
%% separate T1 and T2
fH = [];
fH(1) = figure;
set(gcf,'Position',paau6FigPos)
for iT = 1:2
    colors = get(gca,'ColorOrder');
    subplot(3,2,iT)
    hold on
    for iPAAU = 1:4
        p1 = plot(twin(1):twin(end), groupMean.PAAUT(:,iPAAU,iT));
        set(p1, 'Color', colors(iPAAU,:), 'LineWidth', 4)
    end
    if iT==2
        legend('P-att','P-unatt','A-att','A-unatt')
    end
    for iPAAU = 1:4
        shadedErrorBar(twin(1):twin(end), groupMean.PAAUT(:,iPAAU,iT), groupSte.PAAUT(:,iPAAU,iT), {'color',colors(iPAAU,:),'LineWidth',4}, 1)
    end
    vline(0,'k');
    xlabel('time (ms)')
    ylabel('wavelet amp')
    title(sprintf('T%d',iT))
    
    colors = [trigBlue; trigRed];
    subplot(3,2,iT+2)
    hold on
    for iAU = 1:2
        p1 = plot(twin(1):twin(end), groupMean.AUT(:,iAU,iT));
        set(p1, 'Color', colors(iAU,:), 'LineWidth', 4)
    end
    if iT==2
        legend('att','unatt')
    end
    for iAU = 1:2
        shadedErrorBar(twin(1):twin(end), groupMean.AUT(:,iAU,iT), groupSte.AUT(:,iAU,iT), {'color',colors(iAU,:),'LineWidth',4}, 1)
    end
    vline(0,'k');
    xlabel('time (ms)')
    ylabel('wavelet amp')
    title(sprintf('T%d',iT))
    
    colors = trigColorsPA4([1 4],:);
    subplot(3,2,iT+4)
    hold on
    for iPA = 1:2
        p1 = plot(twin(1):twin(end), groupMean.PAT(:,iPA,iT));
        set(p1, 'Color', colors(iPA,:), 'LineWidth', 4)
    end
    if iT==2
        legend('present','absent')
    end
    for iPA = 1:2
        shadedErrorBar(twin(1):twin(end), groupMean.PAT(:,iPA,iT), groupSte.PAT(:,iPA,iT), {'color',colors(iPA,:),'LineWidth',4}, 1)
    end
    vline(0,'k');
    xlabel('time (ms)')
    ylabel('wavelet amp')
    title(sprintf('T%d',iT))
end
rd_supertitle2(figTitle)

%% combined across T1 and T2
fH(2) = figure;
set(gcf,'Position',paau3FigPos)
colors = get(gca,'ColorOrder');
subplot(3,1,1)
hold on
for iPAAU = 1:4
    p1 = plot(twin(1):twin(end), groupMean.PAAU(:,iPAAU));
    set(p1, 'Color', colors(iPAAU,:), 'LineWidth', 4)
end
legend('P-att','P-unatt','A-att','A-unatt')
for iPAAU = 1:4
    shadedErrorBar(twin(1):twin(end), groupMean.PAAU(:,iPAAU), groupSte.PAAU(:,iPAAU), {'color',colors(iPAAU,:),'LineWidth',4}, 1)
end
vline(0,'k');
xlabel('time (ms)')
ylabel('wavelet amp')
title('T1 & T2')

colors = [trigBlue; trigRed];
subplot(3,1,2)
hold on
for iAU = 1:2
    p1 = plot(twin(1):twin(end), groupMean.AU(:,iAU));
    set(p1, 'Color', colors(iAU,:), 'LineWidth', 4)
end
legend('att','unatt')
for iAU = 1:2
    shadedErrorBar(twin(1):twin(end), groupMean.AU(:,iAU), groupSte.AU(:,iAU), {'color',colors(iAU,:),'LineWidth',4}, 1)
end
vline(0,'k');
xlabel('time (ms)')
ylabel('wavelet amp')
title('T1 & T2')

colors = trigColorsPA4([1 4],:);
subplot(3,1,3)
hold on
for iPA = 1:2
    p1 = plot(twin(1):twin(end), groupMean.PA(:,iPA));
    set(p1, 'Color', colors(iPA,:), 'LineWidth', 4)
end
legend('present','absent')
for iPA = 1:2
    shadedErrorBar(twin(1):twin(end), groupMean.PA(:,iPA), groupSte.PA(:,iPA), {'color',colors(iPA,:),'LineWidth',4}, 1)
end
vline(0,'k');
xlabel('time (ms)')
ylabel('wavelet amp')
title('T1 & T2')

rd_supertitle2(figTitle)

%% all combined
fH(3) = figure;
set(gcf,'Position',tsFigPos)
shadedErrorBar(t, groupMean.ampsAll, groupSte.ampsAll, {'color','k','LineWidth',4})
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
legend('all trials')
xlabel('time (ms)')
ylabel('wavelet amp')
title(figTitle)

%% save
if saveFigs
    figPrefix = sprintf('%s_plot', figStr);
    figNames = {sprintf('%sAmpsPAAU', measure), ...
        sprintf('%sAmpsPAAUT1T2Comb', measure), ...
        sprintf('%sAmpsAll', measure)};
    rd_saveAllFigs(fH, figNames, figPrefix, figDir);
end

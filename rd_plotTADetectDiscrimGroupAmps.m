function rd_plotTADetectDiscrimGroupAmps(A, measure, subjects, groupData, groupMean, groupSte, saveFigs, figDir, figStr)

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

tsFigPos0 = [0 500 850 375];
tf9FigPos = [0 250 1280 580];

eventTimes = A.eventTimes;
trigNames = A.trigNames;
attNames = A.attNames;
paNames = A.PANames;
t = A.t;

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
    case 'w'
        ylims = [0 400];
        diffYLims = [-50 50];
        diffYLimsGroup = [-20 20];
    case 'h'
        ylims = [-5 30];
        diffYLims = [-5 5];
        diffYLimsGroup = [-2.5 2.5];
    otherwise
        error('measure not recognized')
end
xlims = [t(1) t(end)];
ncols = ceil(sqrt(nSubjects));
nrows = ceil(nSubjects/ncols);

allColors.amps = trigColors;
allColors.ampsAtt = trigColorsAtt2;
allColors.ampsPA = trigColorsPA4;

fieldNames = fieldnames(groupData);
nFields = numel(fieldNames);

for iF = 1:nFields
    fieldName = fieldNames{iF};
    figNamesIndiv{iF} = sprintf('%s%s%sIndiv', measure, upper(fieldName(1)), fieldName(2:end));
    figNamesGroup{iF} = sprintf('%s%s%sGroup', measure, upper(fieldName(1)), fieldName(2:end));
end

%% calculate attT2-T1
valsDiff = diff(groupData.ampsAtt);
valsDiffMean = squeeze(mean(valsDiff,3));
valsDiffSte = squeeze(std(valsDiff,0,3)./sqrt(nSubjects));

valsDiffAbsMean = squeeze(mean(abs(valsDiff),3));
valsDiffAbsSte = squeeze(std(abs(valsDiff),0,3)./sqrt(nSubjects));

%% calculate pres-abs ste
t1PA = cat(1, mean(groupData.ampsPA([1 3],:,:)), mean(groupData.ampsPA([2 4],:,:)));
t2PA = cat(1, mean(groupData.ampsPA([1 2],:,:)), mean(groupData.ampsPA([3 4],:,:)));
t1PADiffSte = squeeze(std(diff(t1PA),0,3)./sqrt(nSubjects)); % T1
t2PADiffSte = squeeze(std(diff(t2PA),0,3)./sqrt(nSubjects)); % T2

%% indiv subjects ts
fH = [];
for iF = 1:nFields
    fieldName = fieldNames{iF};
    
    vals = groupData.(fieldName);
    valsMean = groupMean.(fieldName);
    colors = allColors.(fieldName);
    
    fH(iF) = figure;
    set(gcf,'Position',tf9FigPos)
    for iSubject = 1:nSubjects
        subplot(nrows,ncols,iSubject)
        hold on
        if strcmp(fieldName, 'amps')
            for iCond=1:size(valsMean,2)
                plot(t, vals(:,plotOrder(iCond),iSubject), 'color', colors(iCond,:))
            end
        else
            for iCond=1:size(valsMean,1)
                plot(t, vals(iCond,:,iSubject), 'color', colors(iCond,:))
            end
        end
        plot(xlims, [0 0], 'k')
        xlim(xlims)
        ylim(ylims)
        for iEv = 1:numel(eventTimes)
            vline(eventTimes(iEv),'color','k','LineStyle',':');
        end
        if iSubject==1
            xlabel('time (ms)')
            ylabel('amplitude')
        end
        title(und2space(subjects{iSubject}))
    end
    rd_supertitle(figTitle);
    rd_raiseAxis(gca);
end

if saveFigs
    figPrefix = sprintf('%s_plot', figStr);
    rd_saveAllFigs(fH, figNamesIndiv, figPrefix, figDir);
end

%% indiv attT2-attT1
figure
set(gcf,'Position',tf9FigPos)
for iSubject = 1:nSubjects
    subplot(nrows,ncols,iSubject)
    hold on
    plot(t, valsDiff(1,:,iSubject), 'k', 'LineWidth', 2)
    plot(xlims, [0 0], 'k')
    xlim(xlims)
    ylim(diffYLims)
    for iEv = 1:numel(eventTimes)
        vline(eventTimes(iEv),'color','k','LineStyle',':');
    end
    if iSubject==1
        xlabel('time (ms)')
        ylabel('amplitude difference (T2-T1)')
    end
    title(und2space(subjects{iSubject}))
end
rd_supertitle(figTitle);
rd_raiseAxis(gca);

if saveFigs
    figPrefix = sprintf('%s_plot', figStr);
    figName = sprintf('%sAmpsAttDiffIndiv', measure);
    rd_saveAllFigs(gcf, {figName}, figPrefix, figDir);
end

%% group
fH = [];
for iF = 1:nFields
    fieldName = fieldNames{iF};
    
    valsMean = groupMean.(fieldName);
    valsSte = groupSte.(fieldName);
    colors = allColors.(fieldName);
    
    if strcmp(fieldName, 'amps')
        valsMean = valsMean(:,plotOrder)';
        valsSte = valsSte(:,plotOrder)';
    end
    
    fH(iF) = figure;
    set(gcf,'Position',tsFigPos0);
    hold on
    for iCond=1:size(valsMean,1)
        shadedErrorBar(t, valsMean(iCond,:), valsSte(iCond,:), {'color', colors(iCond,:), 'LineWidth', 3}, 1)
    end
%     plot(xlims, [0 0], 'k')
    for iEv = 1:numel(eventTimes)
        vline(eventTimes(iEv),'color','k','LineStyle',':');
    end
    xlim(xlims)
    xlabel('time (ms)')
    ylabel('amplitude')
    title(figTitle)
end

if saveFigs
    figPrefix = sprintf('%s_plot', figStr);
    rd_saveAllFigs(fH, figNamesGroup, figPrefix, figDir);
end

%% group attT1 vs attT2, separated by PA
valsMean = groupMean.amps';
valsSte = groupSte.amps';
colors = allColors.ampsAtt;
conds = 1:2:nTrigs-1;

fH = [];
fH(1) = figure;
for iCond = 1:numel(conds)
    subplot(2,2,iCond)
    hold on
    shadedErrorBar(t, valsMean(conds(iCond),:), valsSte(conds(iCond),:), {'color', colors(1,:), 'LineWidth', 3}, 1)
    shadedErrorBar(t, valsMean(conds(iCond)+1,:), valsSte(conds(iCond)+1,:), {'color', colors(2,:), 'LineWidth', 3}, 1)
    for iEv = 1:numel(eventTimes)
        vline(eventTimes(iEv),'color','k','LineStyle',':');
    end
    xlim(xlims)
    xlabel('time (ms)')
    ylabel('amplitude')
    title(paNames{iCond})
end

if saveFigs
    figPrefix = sprintf('%s_plot', figStr);
    figNames = {sprintf('%sAmpsAttByPA', measure)};
    rd_saveAllFigs(fH, figNames, figPrefix, figDir);
end

%% group attT2-attT1 with ste error bars
valsMean = groupMean.ampsAtt;
colors = allColors.ampsAtt;

fH = [];
fH(1) = figure;
set(gcf,'Position',tsFigPos0);
hold on
for iCond=1:2
    shadedErrorBar(t, valsMean(iCond,:), valsDiffSte, {'color', colors(iCond,:), 'LineWidth', 3}, 1)
end
% plot(t, groupMean.amps(:,end), 'color', allColors.amps(end,:), 'LineWidth', 1.5)
% plot(xlims, [0 0], 'k')
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'color','k','LineStyle',':');
end
xlim(xlims)
xlabel('time (ms)')
ylabel('amplitude')
title(figTitle)

fH(2) = figure;
set(gcf,'Position',tsFigPos0);
hold on
shadedErrorBar(t, valsDiffMean, valsDiffSte, 'k', 1)
plot(xlims, [0 0], 'k')
ylim(diffYLimsGroup)
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'color','k','LineStyle',':');
end
xlim(xlims)
xlabel('time (ms)')
ylabel('amplitude difference (T2-T1)')
title(figTitle)

fH(3) = figure;
set(gcf,'Position',tsFigPos0);
hold on
shadedErrorBar(t, valsDiffAbsMean, valsDiffAbsSte, 'k', 1)
plot(xlims, [0 0], 'k')
ylim(diffYLimsGroup)
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'color','k','LineStyle',':');
end
xlim(xlims)
xlabel('time (ms)')
ylabel('|amplitude difference (T2-T1)|')
title(figTitle)

if saveFigs
    figPrefix = sprintf('%s_plot', figStr);
    figNames = {sprintf('%sAmpsAttGroupSte', measure), ...
        sprintf('%sAmpsAttDiffGroupSte', measure), ...
        sprintf('%sAmpsAttDiffAbsGroupSte', measure)};
    rd_saveAllFigs(fH, figNames, figPrefix, figDir);
end

%% group pres-abs with ste error bars
colors = allColors.ampsPA([1 4],:);
ylims = [135 235];
twin = [-200 700];
t1Tidx = find(t==eventTimes(3)+twin(1)):find(t==eventTimes(3)+twin(2));
t2Tidx = find(t==eventTimes(4)+twin(1)):find(t==eventTimes(4)+twin(2));

fH = [];
fH(1) = figure;
set(gcf,'Position',tsFigPos0);
subplot(1,2,1)
hold on
for iPA=1:2
    shadedErrorBar(twin(1):twin(end), mean(t1PA(iPA,t1Tidx,:),3), t1PADiffSte(t1Tidx), {'color', colors(iPA,:), 'LineWidth', 3}, 1)
%     shadedErrorBar(t(t1Tidx), mean(t1PA(iPA,t1Tidx,:),3), t1PADiffSte(t1Tidx), {'color', colors(iPA,:), 'LineWidth', 3}, 1)
end
ylim(ylims)
vline(0,'color','k','LineStyle',':'); %%%%
% for iEv = 1:numel(eventTimes)
%     vline(eventTimes(iEv),'color','k','LineStyle',':');
% end
xlim(twin)
% xlim([t(t1Tidx(1)) t(t1Tidx(end))])
xlabel('time (ms)')
ylabel('amplitude')
title('T1')

subplot(1,2,2)
hold on
for iPA=1:2
    shadedErrorBar(twin(1):twin(end), mean(t2PA(iPA,t2Tidx,:),3), t2PADiffSte(t2Tidx), {'color', colors(iPA,:), 'LineWidth', 3}, 1)
%     shadedErrorBar(t(t2Tidx), mean(t2PA(iPA,t2Tidx,:),3), t2PADiffSte(t2Tidx), {'color', colors(iPA,:), 'LineWidth', 3}, 1)
end
ylim(ylims)
vline(0,'color','k','LineStyle',':'); %%%%
% for iEv = 1:numel(eventTimes)
%     vline(eventTimes(iEv),'color','k','LineStyle',':');
% end
xlim(twin)
% xlim([t(t2Tidx(1)) t(t2Tidx(end))])
xlabel('time (ms)')
ylabel('amplitude')
title('T2')

rd_supertitle(figTitle)

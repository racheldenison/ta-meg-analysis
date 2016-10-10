function rd_plotTADetectDiscrimGroupAmps2(A, measure, subjects, groupData, groupMean, groupSte, saveFigs, figDir, figStr)

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
tsFigPos0 = [0 500 850 375];
tf9FigPos = [0 250 1280 580];
paau3FigPos = [800 30 450 830];
paau6FigPos = [0 90 1000 800];

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
        switch A.normalizeOption
            case 'none'
                ylims = [0 400];
                diffYLims = [-50 50];
                diffYLimsGroup = [-20 20];
            case 'commonBaseline'
                ylims = [-1 2];
                diffYLims = [-.8 .8];
                diffYLimsGroup = [-.3 .3];               
            otherwise
                ylims = [0 2];
                diffYLims = [-.8 .8];
                diffYLimsGroup = [-.3 .3];
        end
    case 'h'
        ylims = [-5 30];
        diffYLims = [-5 5];
        diffYLimsGroup = [-2.5 2.5];
    case 'w-single'
        switch A.normalizeOption
            case 'none'
                ylims = [300 700];
                diffYLims = [-1.5 1.5];
                diffYLimsGroup = [-1 1];
            otherwise
                ylims = [.6 1.4];
                diffYLims = [-.8 .8];
                diffYLimsGroup = [-.3 .3];
        end
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

%% permutation test: any difference between att T1 and att T2?
% shuffle condition labels to generate null distribution of attDiff
shuffle = 0;
if shuffle
    nShuffles = 1000;
    vals = groupData.ampsAtt;
    
    for iShuffle = 1:nShuffles
        for iT = 1:2
            attLabel = randi(2, 1, nSubjects); % att T1, att T2
            for iS = 1:nSubjects
                shuffleData(:,1,iS,iShuffle) = vals(attLabel(iS),:,iS);
                shuffleData(:,2,iS,iShuffle) = vals(3-attLabel(iS),:,iS);
            end
        end
    end
    shuffleAttDiff = squeeze(diff(shuffleData,1,2)); % [time subject shuffle]
    shuffleAttDiffMean = squeeze(mean(shuffleAttDiff,2)); % [time shuffle]
    
    ci = prctile(shuffleAttDiffMean,[2.5 97.5],2);
    
    figure
    hold on
    plot(t,valsDiffMean)
    plot(t,valsDiffCI,'g') % time series shuffle
    plot(t,ci) % condition label flip
end

%% indiv subjects ts
fH = [];
for iF = 1:nFields
    fieldName = fieldNames{iF};
    
    vals = groupData.(fieldName);
    valsMean = groupMean.(fieldName);
    colors = allColors.(fieldName);
%     if isfield(allColors, fieldName)
%         colors = allColors.(fieldName);
%     else
%         break % do not plot
%     end
    
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
%     if isfield(allColors, fieldName)
%         colors = allColors.(fieldName);
%     else
%         break % do not plot
%     end
    
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


%% PAAU
%% calculate PAAU
twin = [-600 600];
t1Tidx = find(t==eventTimes(3)+twin(1)):find(t==eventTimes(3)+twin(2));
t2Tidx = find(t==eventTimes(4)+twin(1)):find(t==eventTimes(4)+twin(2));
twindow = twin(1):twin(end);

% calculate pres/abs x att/unattend for each target, groupData
groupData.PAAUT(:,1,1,:) = mean(groupData.amps(t1Tidx,[1 5],:),2); % present/attended
groupData.PAAUT(:,2,1,:) = mean(groupData.amps(t1Tidx,[2 6],:),2); % present/unattended
groupData.PAAUT(:,3,1,:) = mean(groupData.amps(t1Tidx,[3 7],:),2); % absent/attended
groupData.PAAUT(:,4,1,:) = mean(groupData.amps(t1Tidx,[4 8],:),2); % absent/unattended

groupData.PAAUT(:,1,2,:) = mean(groupData.amps(t2Tidx,[2 4],:),2);
groupData.PAAUT(:,2,2,:) = mean(groupData.amps(t2Tidx,[1 3],:),2);
groupData.PAAUT(:,3,2,:) = mean(groupData.amps(t2Tidx,[6 8],:),2);
groupData.PAAUT(:,4,2,:) = mean(groupData.amps(t2Tidx,[5 7],:),2);

% present vs. absent and attended vs. unattended
for iT = 1:2
    groupData.PAT(:,1,iT,:) = mean(groupData.PAAUT(:,[1 2],iT,:),2); % present
    groupData.PAT(:,2,iT,:) = mean(groupData.PAAUT(:,[3 4],iT,:),2); % absent

    groupData.AUT(:,1,iT,:) = mean(groupData.PAAUT(:,[1 3],iT,:),2); % attended
    groupData.AUT(:,2,iT,:) = mean(groupData.PAAUT(:,[2 4],iT,:),2); % unattended 
end

% combining across T1 and T2
for iPAAU = 1:4
    groupData.PAAU(:,iPAAU,:) = mean(groupData.PAAUT(:,iPAAU,[1 2],:),3);
end

groupData.PA(:,1,:) = mean(groupData.PAAU(:,[1 2],:),2); % present
groupData.PA(:,2,:) = mean(groupData.PAAU(:,[3 4],:),2); % absent

groupData.AU(:,1,:) = mean(groupData.PAAU(:,[1 3],:),2); % attended
groupData.AU(:,2,:) = mean(groupData.PAAU(:,[2 4],:),2); % unattended

% finally, ampsAll
groupData.ampsAll = squeeze(mean(groupData.amps(:,1:end-1,:),2));

% group means and ste
fieldNames = {'PAAUT','PAT','AUT','PAAU','PA','AU','ampsAll'};
nFields = numel(fieldNames);
for iF = 1:nFields
    fieldName = fieldNames{iF};
    vals = groupData.(fieldName);
    sdim = numel(size(vals)); % subject dimension
    groupMean.(fieldName) = mean(vals, sdim);
    groupSte.(fieldName) = std(vals, 0, sdim)./sqrt(nSubjects);
end

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

%% enhancers/suppressers
% calculate mean amp over a window
twinm = [-250 50]; % [-200 200]
tidx = find(twindow==twinm(1)):find(twindow==twinm(2));
winamp(:,1,:) = mean(groupData.PAAUT(tidx,:,1,:),1); % T1
winamp(:,2,:) = mean(groupData.PAAUT(tidx,:,2,:),1); % T2
winampAU(1,:) = mean(winamp(3,:,:),2); % absent/attended
winampAU(2,:) = mean(winamp(4,:,:),2); % absent/unattended
enhancers = subjects(diff(winampAU)<0);
suppressers = subjects(diff(winampAU)>0);

% winampAU2(1,:,:) = mean(winamp([1 3],:,:));
% winampAU2(2,:,:) = mean(winamp([2 4],:,:));
% winampAU2GroupMean = mean(winampAU2,3);
% winampAU2GroupSte = std(winampAU2,0,3)/sqrt(nSubjects);
% figure
% barweb(winampAU2GroupMean',winampAU2GroupSte');
% winampAU2N = normalizeDC(winampAU2);
% wnmean = mean(winampAU2N,3);
% wnste = std(winampAU2N,0,3)./sqrt(nSubjects);
% figure
% barweb(wnmean',wnste');
% 
% % effect size
% winampAU2Diff = winampAU2(1,:,:)-winampAU2(2,:,:);
% cohenD = mean(winampAU2Diff,3)./std(winampAU2Diff,0,3);

% att-unatt consistency (correlation across subjects for P and A)
a1(:,1,:,:) = squeeze(groupData.PAAUT(:,1,:,:)-groupData.PAAUT(:,2,:,:)); % time x P/A x target x subject
a1(:,2,:,:) = squeeze(groupData.PAAUT(:,3,:,:)-groupData.PAAUT(:,4,:,:));

for i = 1:numel(twindow)
    temp = corr(squeeze(a1(i,:,1,:))'); % T1
    a1corr(i) = temp(1,2);
    temp = corr(squeeze(a1(i,:,2,:))'); % T2
    a2corr(i) = temp(1,2);
end

figure
hold on
plot(twindow,a1corr)
plot(twindow,a2corr,'r')
plot(twindow([1 end]),[0 0],'k:')
vline(0,'color','k','LineStyle',':');
xlabel('time (ms)')
ylabel('correlation between present att-unatt and absent att-unatt')
legend('T1','T2')

% figure
% for i = 1:16
%     subplot(4,4,i)
%     bar(winamp(:,:,i)')
% end
% legend('P-att','P-unatt','A-att','A-unatt')
% 
% figure
% bar(winampAU')
% legend('A-att','A-unatt')

if saveFigs
    figPrefix = sprintf('%s_plot', figStr);
    rd_saveAllFigs(gcf, {'wAmpsPresVsAbsAUDiffCorr'}, figPrefix, figDir);
end


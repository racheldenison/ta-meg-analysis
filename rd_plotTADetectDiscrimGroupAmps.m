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
    case 'w-single'
        ylims = [300 700];
        diffYLims = [-1.5 1.5];
        diffYLimsGroup = [-1 1];
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

% generate null distribution by shuffling each subject's time series in time
nSamples = 1000;
for iSample = 1:nSamples
    for iSubject=1:nSubjects
        valsDiffShuffled(:,iSubject,iSample) = squeeze(shuffle(valsDiff(:,:,iSubject),2));
    end
end
valsDiffCI = prctile(squeeze(nanmean(valsDiffShuffled,2)),[2.5 97.5],2);

% absolute value
valsDiffShuffledAbsMean = squeeze(nanmean(abs(valsDiffShuffled),2));
% valsDiffAbsCI = prctile(valsDiffShuffledAbsMean,[2.5 97.5],2);
valsDiffAbsCI = prctile(valsDiffShuffledAbsMean,95,2);

%% permutation test: any difference between att T1 and att T2?
% shuffle condition labels to generate null distribution of attDiff
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

% figure
% hold on
% plot(t,valsDiffMean)
% plot(t,valsDiffCI,'g') % time series shuffle
% plot(t,ci) % condition label flip

%% calculate pres-abs ste
t1PA = cat(1, mean(groupData.ampsPA([1 3],:,:)), mean(groupData.ampsPA([2 4],:,:)));
t2PA = cat(1, mean(groupData.ampsPA([1 2],:,:)), mean(groupData.ampsPA([3 4],:,:)));
t1PADiffSte = squeeze(std(diff(t1PA),0,3)./sqrt(nSubjects)); % T1
t2PADiffSte = squeeze(std(diff(t2PA),0,3)./sqrt(nSubjects)); % T2

%% calculate pres/abs x att/unattend for each target, groupMean
t1PAAU(:,1) = mean(groupMean.amps(:,[1 5]),2); % present/attended
t1PAAU(:,2) = mean(groupMean.amps(:,[2 6]),2); % present/unattended
t1PAAU(:,3) = mean(groupMean.amps(:,[3 7]),2); % absent/attended
t1PAAU(:,4) = mean(groupMean.amps(:,[4 8]),2); % absent/unattended

t2PAAU(:,1) = mean(groupMean.amps(:,[2 4]),2);
t2PAAU(:,2) = mean(groupMean.amps(:,[1 3]),2);
t2PAAU(:,3) = mean(groupMean.amps(:,[6 8]),2);
t2PAAU(:,4) = mean(groupMean.amps(:,[5 7]),2);

%% calculate pres/abs x att/unattend for each target, groupData
t1PAAUData(:,1,:) = mean(groupData.amps(:,[1 5],:),2); % present/attended
t1PAAUData(:,2,:) = mean(groupData.amps(:,[2 6],:),2); % present/unattended
t1PAAUData(:,3,:) = mean(groupData.amps(:,[3 7],:),2); % absent/attended
t1PAAUData(:,4,:) = mean(groupData.amps(:,[4 8],:),2); % absent/unattended

t2PAAUData(:,1,:) = mean(groupData.amps(:,[2 4],:),2);
t2PAAUData(:,2,:) = mean(groupData.amps(:,[1 3],:),2);
t2PAAUData(:,3,:) = mean(groupData.amps(:,[6 8],:),2);
t2PAAUData(:,4,:) = mean(groupData.amps(:,[5 7],:),2);

t1PAAUSte = std(t1PAAUData,0,3)/sqrt(nSubjects); 
t2PAAUSte = std(t2PAAUData,0,3)/sqrt(nSubjects);

%% calculate mean amp over a window
twin = [-100 100]; % [-200 200]
t1Tidx = find(t==eventTimes(3)+twin(1)):find(t==eventTimes(3)+twin(2));
t2Tidx = find(t==eventTimes(4)+twin(1)):find(t==eventTimes(4)+twin(2));
winamp(:,1,:) = mean(t1PAAUData(t1Tidx,:,:),1); % T1
winamp(:,2,:) = mean(t2PAAUData(t2Tidx,:,:),1); % T2
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
a1(:,:,1) = t1PAAUData(:,1,:)-t1PAAUData(:,2,:);
a1(:,:,2) = t1PAAUData(:,3,:)-t1PAAUData(:,4,:);

a2(:,:,1) = t2PAAUData(:,1,:)-t2PAAUData(:,2,:);
a2(:,:,2) = t2PAAUData(:,3,:)-t2PAAUData(:,4,:);

for i = 1:numel(t)
    temp = corr(squeeze(a1(i,:,:)));
    a1corr(i) = temp(1,2);
    temp = corr(squeeze(a2(i,:,:)));
    a2corr(i) = temp(1,2);
end

twin = [-200 700];
t1Tidx = find(t==eventTimes(3)+twin(1)):find(t==eventTimes(3)+twin(2));
t2Tidx = find(t==eventTimes(4)+twin(1)):find(t==eventTimes(4)+twin(2));

a1corr = a1corr(t1Tidx);
a2corr = a2corr(t2Tidx);

figure
hold on
plot(twin(1):twin(end),a1corr)
plot(twin(1):twin(end),a2corr,'r')
plot(twin,[0 0],'k:')
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

% PA
fH = [];
fH(1) = figure;
set(gcf,'Position',tsFigPos0);
subplot(1,2,1)
hold on
for iPA=1:2
    shadedErrorBar(twin(1):twin(end), mean(t1PA(iPA,t1Tidx,:),3), t1PADiffSte(t1Tidx), {'color', colors(iPA,:), 'LineWidth', 3}, 1)
end
ylim(ylims)
vline(0,'color','k','LineStyle',':');
xlim(twin)
% xlim([t(t1Tidx(1)) t(t1Tidx(end))])
xlabel('time (ms)')
ylabel('amplitude')
title('T1')

subplot(1,2,2)
hold on
for iPA=1:2
    shadedErrorBar(twin(1):twin(end), mean(t2PA(iPA,t2Tidx,:),3), t2PADiffSte(t2Tidx), {'color', colors(iPA,:), 'LineWidth', 3}, 1)
end
ylim(ylims)
vline(0,'color','k','LineStyle',':');
xlim(twin)
% xlim([t(t2Tidx(1)) t(t2Tidx(end))])
xlabel('time (ms)')
ylabel('amplitude')
title('T2')

rd_supertitle(figTitle)

% PAAU
fH(2) = figure;
colors = get(gca,'ColorOrder');
set(gcf,'Position',tsFigPos0);
subplot(1,2,1)
hold on
plot(twin(1):twin(end), t1PAAU(t1Tidx,:))
for iPAAU = 1:4
    shadedErrorBar(twin(1):twin(end), t1PAAU(t1Tidx,iPAAU), t1PAAUSte(t1Tidx,iPAAU), {'color', colors(iPAAU,:), 'LineWidth', 3}, 1)
end
% ylim(ylims)
vline(0,'color','k','LineStyle',':'); 
xlim(twin)
% xlim([t(t1Tidx(1)) t(t1Tidx(end))])
xlabel('time (ms)')
ylabel('amplitude')
title('T1')

subplot(1,2,2)
hold on
plot(twin(1):twin(end), t2PAAU(t2Tidx,:))
for iPAAU = 1:4
    shadedErrorBar(twin(1):twin(end), t2PAAU(t2Tidx,iPAAU), t2PAAUSte(t2Tidx,iPAAU), {'color', colors(iPAAU,:), 'LineWidth', 3}, 1)
end
% ylim(ylims)
vline(0,'color','k','LineStyle',':');
xlim(twin)
% xlim([t(t2Tidx(1)) t(t2Tidx(end))])
xlabel('time (ms)')
ylabel('amplitude')
title('T2')
legend('P-att','P-unatt','A-att','A-unatt')

rd_supertitle(figTitle)

paau = (t1PAAU(t1Tidx,:) + t2PAAU(t2Tidx,:))/2;
fH(3) = figure;
plot(twin(1):twin(end),paau)
vline(0,'color','k','LineStyle',':');
legend('P-att','P-unatt','A-att','A-unatt')
xlabel('time (ms)')
ylabel('amplitude')
title('T1 & T2')
box off

paauData = (t1PAAUData(t1Tidx,:,:) + t2PAAUData(t2Tidx,:,:))/2;
paauPresAUDiff = squeeze(paauData(:,1,:) - paauData(:,2,:));
fH(4) = figure;
set(gcf,'Position',tf9FigPos)
for iSubject = 1:nSubjects
    subplot(nrows,ncols,iSubject)
    hold on
    plot(twin(1):twin(end), paauData(:,1:2,iSubject), 'LineWidth', 2)
%     xlim(xlims)
%     ylim(diffYLims)
    vline(0,'color','k','LineStyle',':');
    if iSubject==1
        xlabel('time (ms)')
%         ylabel('amplitude difference (T2-T1)')
    end
    title(und2space(subjects{iSubject}))
end
legend('P-att','P-unatt')

fH(5) = figure;
hold on
plot(twin([1 end]), [0 0], 'k:')
shadedErrorBar(twin(1):twin(end), mean(paauPresAUDiff,2), std(paauPresAUDiff,0,2)/sqrt(nSubjects), {'color', 'k', 'LineWidth', 3}, 1)
vline(0,'color','k','LineStyle',':');
xlabel('time (ms)')
ylabel('amplitude difference (att-unatt)')
title('target present trials')

if saveFigs
    figPrefix = sprintf('%s_plot', figStr);
    figNames = {sprintf('%sAmpsPAByTargetGroup', measure), ...
        sprintf('%sAmpsPAAUGroupSte', measure), ...
        sprintf('%sAmpsPAAUT1T2AveGroup', measure), ...
        sprintf('%sAmpsPAUIndiv', measure), ...
        sprintf('%sAmpsPAUDiffGroup', measure)};
    rd_saveAllFigs(fH, figNames, figPrefix, figDir);
end

%% paau separately for T1 and T2
paauDataT1T2(:,:,:,1) = t1PAAUData(t1Tidx,:,:); 
paauDataT1T2(:,:,:,2) = t2PAAUData(t2Tidx,:,:);
paauPresAUDiffT1T2 = squeeze(paauDataT1T2(:,1,:,:) - paauDataT1T2(:,2,:,:));

for iT = 1:2
    fH(5+iT) = figure;
    set(gcf,'Position',tf9FigPos)
    for iSubject = 1:nSubjects
        subplot(nrows,ncols,iSubject)
        hold on
        plot(twin(1):twin(end), paauDataT1T2(:,1:2,iSubject,iT), 'LineWidth', 2)
        %     xlim(xlims)
        %     ylim(diffYLims)
        vline(0,'color','k','LineStyle',':');
        if iSubject==1
            xlabel('time (ms)')
            %         ylabel('amplitude difference (T2-T1)')
        end
        title(und2space(subjects{iSubject}))
    end
    legend('P-att','P-unatt')
    rd_supertitle2(sprintf('T%d', iT))
end

fH(8) = figure;
for iT = 1:2
    subplot(1,2,iT)
    hold on
    plot(twin([1 end]), [0 0], 'k:')
    shadedErrorBar(twin(1):twin(end), mean(paauPresAUDiffT1T2(:,:,iT),2), std(paauPresAUDiffT1T2(:,:,iT),0,2)/sqrt(nSubjects), {'color', 'k', 'LineWidth', 3}, 1)
    vline(0,'color','k','LineStyle',':');
    xlabel('time (ms)')
    ylabel('amplitude difference (att-unatt)')
    title(sprintf('T%d, target present trials', iT))
end


% rd_plotTADetectDiscrimGroupHMFACR.m

%% setup
% trialSelections = {'detectHit','detectMiss','detectFA','detectCR'};
trialSelections = {'correct','incorrect'};
respTarget = 'T1Resp';

%% get data
for i = 1:numel(trialSelections)
    selectionStr = sprintf('topChannels5_%sTrials%s', trialSelections{i}, respTarget);
    [groupDataTS(i), groupMeanTS(i), groupSteTS(i), A] = rd_plotTADetectDiscrimGroup('ts-single', selectionStr, 'none');
end

for i = 1:numel(trialSelections)
    selectionStr = sprintf('topChannels5_%sTrials%s', trialSelections{i}, respTarget);
    [groupDataW(i), groupMeanW(i), groupSteW(i), A]  = rd_plotTADetectDiscrimGroup('w-single', selectionStr, 'stim');
end

%% organize hmfc data
switch respTarget
    case 'T1Resp'
        iT = 1;
    case 'T2Resp'
        iT = 2;
    otherwise
        error('respTarget not recognized')
end
for i = 1:numel(trialSelections)
    hmfcMean.nTrialsPerCond(:,i) = groupMeanTS(i).nTrialsPerCond;
    hmfcSte.nTrialsPerCond(:,i) = groupSteTS(i).nTrialsPerCond;
    
    hmfcMean.PAAUT(:,:,i) = groupMeanW(i).PAAUT(:,:,iT);
    hmfcSte.PAAUT(:,:,i) = groupSteW(i).PAAUT(:,:,iT);
    
    hmfcMean.PAT(:,:,i) = groupMeanW(i).PAT(:,:,iT);
    hmfcSte.PAT(:,:,i) = groupSteW(i).PAT(:,:,iT);
end

%% plot figs
twin = A.wtwin;
colors = get(gca,'ColorOrder');
lw = [4 2 4 2];

%% nTrials
figure
hold on
bar(hmfcMean.nTrialsPerCond)
plot([1 numel(A.trigNames)],[56 56],'--k')
ylim([0 60])
ylabel('number of trials')
set(gca,'XTick', 1:numel(A.trigNames))
set(gca,'XTickLabel',A.trigNames)
rotateXLabels(gca,45)
legend(trialSelections)
colormap(colors(1:4,:))

%% PAAUT
ylims = [.9 1.15];
figure
subplot(1,2,1)
hold on
plot(twin(1):twin(end), squeeze(hmfcMean.PAAUT(:,1,:)),'LineWidth',4)
plot(twin(1):twin(end), squeeze(hmfcMean.PAAUT(:,2,:)),'LineWidth',2)
for i = 1:2
    for j = 1:2
        shadedErrorBar(twin(1):twin(end), squeeze(hmfcMean.PAAUT(:,i,j)), ...
            squeeze(hmfcSte.PAAUT(:,i,j)), {'color',colors(j,:),'LineWidth',lw(i)}, 1)
    end
end
xlim(twin)
ylim(ylims);
vline(0,'k');
xlabel('time (ms)')
ylabel('wavelet amp')
title('present')
subplot(1,2,2)
hold on
plot(twin(1):twin(end), squeeze(hmfcMean.PAAUT(:,3,:)),'LineWidth',4)
plot(twin(1):twin(end), squeeze(hmfcMean.PAAUT(:,4,:)),'LineWidth',2)
for i = 3:4
    for j = 3:4
        shadedErrorBar(twin(1):twin(end), squeeze(hmfcMean.PAAUT(:,i,j)), ...
            squeeze(hmfcSte.PAAUT(:,i,j)), {'color',colors(j,:),'LineWidth',lw(i)}, 1)
    end
end
xlim(twin)
ylim(ylims);
vline(0,'k');
xlabel('time (ms)')
ylabel('wavelet amp')
title('absent')
legend(trialSelections)

%% PAT
ylims = [.9 1.15];
figure
subplot(1,2,1)
hold on
plot(twin(1):twin(end), squeeze(hmfcMean.PAT(:,1,:)),'LineWidth',4)
for j = 1:2
    shadedErrorBar(twin(1):twin(end), squeeze(hmfcMean.PAAUT(:,1,j)), ...
        squeeze(hmfcSte.PAAUT(:,1,j)), {'color',colors(j,:),'LineWidth',4}, 1)
end
xlim(twin)
ylim(ylims);
vline(0,'k');
xlabel('time (ms)')
ylabel('wavelet amp')
title('present')
subplot(1,2,2)
hold on
plot(twin(1):twin(end), squeeze(hmfcMean.PAT(:,2,:)),'LineWidth',4)
for j = 1:2
    shadedErrorBar(twin(1):twin(end), squeeze(hmfcMean.PAT(:,2,j)), ...
        squeeze(hmfcSte.PAT(:,2,j)), {'color',colors(j,:),'LineWidth',4}, 1)
end
xlim(twin)
ylim(ylims);
vline(0,'k');
xlabel('time (ms)')
ylabel('wavelet amp')
title('absent')
legend(trialSelections)

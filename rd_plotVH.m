% rd_plotVH.m

%%
exptDir = pathToTANoise('MEG');
analysisDir = sprintf('%s/Group/mat', exptDir);

%% load data
measure = 'itpc-single';
selectionStr = 'topChannels5_allTrials';
normalizeOption = 'none';
% also set exptType and ssvefFreq in function

[groupData, groupMean, groupSte, A] = rd_plotTADetectDiscrimGroup(measure, selectionStr, normalizeOption);

%% setup
targetNames = {'T1','T2'};
attNames = {'att','unatt'};

nT = numel(targetNames);
nAtt = numel(attNames);

twin = A.wtwin;
tt = twin(1):twin(end);

nTrigs = size(A.ampsMean,3);
nSubjects = size(groupData.amps,3)/2; % /2 assuming 2 sessions per subject

% colors
extendedMap = flipud(lbmap(nTrigs-1+4,'RedBlue'));
selectedMap = extendedMap([1:(nTrigs-1)/2 (end-(nTrigs-1)/2)+1:end],:);
trigColors = [selectedMap; 0 0 0];
trigBlue = mean(selectedMap(1:(nTrigs-1)/2,:));
trigRed = mean(selectedMap((end-(nTrigs-1)/2)+1:end,:));
trigColorsAtt2 = [trigBlue; trigRed];
trigColorsPA4 = [.52 .37 .75; .31 .74 .40; .27 .51 .84; 1.0 .57 .22];

allColors.amps = trigColors;
allColors.ampsAtt = trigColorsAtt2;
allColors.ampsPA = trigColorsPA4;

%% calculate V-H and att-unatt combinations
vhAve(:,1,:,:) = squeeze(mean(groupData.PAAUT(:,1:2,:,:),2)); % v
vhAve(:,2,:,:) = squeeze(mean(groupData.PAAUT(:,3:4,:,:),2)); % h

vhDiff(:,1,:,:) = groupData.PAAUT(:,1,:,:) - groupData.PAAUT(:,3,:,:);
vhDiff(:,2,:,:) = groupData.PAAUT(:,2,:,:) - groupData.PAAUT(:,4,:,:);

vhAveS = [];
vhDiffS = [];
for iS = 1:nSubjects
    vhAveS(:,:,:,iS) = (vhAve(:,:,:,iS*2-1) + vhAve(:,:,:,iS*2))/2;
    vhDiffS(:,:,:,iS) = (vhDiff(:,:,:,iS*2-1) + vhDiff(:,:,:,iS*2))/2;
end
vhAveSVHDiff = squeeze(vhAveS(:,1,:,:) - vhAveS(:,2,:,:));
vhAveSTAve = squeeze(mean(vhAveS,3));
vhAveSTAveVHDiff = squeeze(vhAveSTAve(:,1,:) - vhAveSTAve(:,2,:));

vhDiffSAttDiff = squeeze(vhDiffS(:,1,:,:) - vhDiffS(:,2,:,:));
vhDiffSTAve = squeeze(mean(vhDiffS,3));
vhDiffSTAveAttDiff = squeeze(vhDiffSTAve(:,1,:) - vhDiffSTAve(:,2,:));

%% V and H overlaid, group error bars
vhColors = allColors.ampsPA([1 4],:);
ylims = [.15 .3];
figure
for iT = 1:nT
    subplot(nT+1,1,iT)
    hold on
    for iVH = 1:2
        m = squeeze(mean(vhAveS(:,iVH,iT,:),4));
%         se = squeeze(std(vhAveS(:,iVH,iT,:),0,4)/sqrt(nSubjects));
        se = squeeze(std(vhAveSVHDiff(:,iVH,:),0,3)/sqrt(nSubjects));
        shadedErrorBar(tt, m, se, ...
            {'Color',vhColors(iVH,:)}, 1);
    end
    ylim(ylims)
    vline(0,'k')
    title(targetNames{iT})
end
subplot(nT+1,1,nT+1)
hold on
for iVH = 1:2
    m = squeeze(mean(vhAveSTAve(:,iVH,:),3));
%     se = squeeze(std(vhAveSTAve(:,iVH,:),0,3)/sqrt(nSubjects));
    se = squeeze(std(vhAveSTAveVHDiff,0,2)/sqrt(nSubjects));
    p1(iVH) = plot(tt, m, 'Color', vhColors(iVH,:));
    shadedErrorBar(tt, m, se, ...
        {'Color',vhColors(iVH,:)}, 1);
end
ylim(ylims)
vline(0,'k')
title('T1 & T2')
xlabel('Time (ms)')
ylabel('ITPC')
legend(p1, {'vertical','horizontal'})


%% V-H
for iT = 1:nT
    for iAtt = 1:nAtt
        figure
        hold on
        plot(tt, squeeze(vhDiffS(:,iAtt,iT,:)))
        plot(tt, squeeze(mean(vhDiffS(:,iAtt,iT,:),4)), '-k', 'LineWidth',4)
        plot(twin,[0 0],'k')
        vline(0,'k')
        
        title(sprintf('%s, V-%s - H-%s',targetNames{iT},attNames{iAtt},attNames{iAtt}))
    end
end

%% V-H, att and unatt overlaid
for iT = 1:nT
    figure
    hold on
    for iAtt = 1:nAtt
        plot(tt, squeeze(vhDiffS(:,iAtt,iT,:)), 'Color',allColors.ampsAtt(iAtt,:))
    end
    for iAtt = 1:nAtt
        p1(iAtt) = plot(tt, squeeze(mean(vhDiffS(:,iAtt,iT,:),4)), 'Color',allColors.ampsAtt(iAtt,:),'LineWidth',4);
    end
    plot(twin,[0 0],'k')
    vline(0,'k')
    legend(p1, attNames)
    title(sprintf('%s, V-H',targetNames{iT}))
end

%% V-H, att and unatt overlaid, group error bars
ylims = [-.065 .065];
figure
for iT = 1:nT
    subplot(nT+1,1,iT)
    hold on
    for iAtt = 1:nAtt
        m = squeeze(mean(vhDiffS(:,iAtt,iT,:),4));
%         se = squeeze(std(vhDiffS(:,iAtt,iT,:),0,4)/sqrt(nSubjects));
        se = squeeze(std(vhDiffSAttDiff(:,iVH,:),0,3)/sqrt(nSubjects));
        shadedErrorBar(tt, m, se, ...
            {'Color',allColors.ampsAtt(iAtt,:)}, 1);
    end
    ylim(ylims)
    vline(0,'k')
    plot(tt([1 end]),[0 0],'k','LineWidth',1)
    title(targetNames{iT})
end
subplot(nT+1,1,nT+1)
hold on
p1 = plot(tt,mean(vhDiffSTAve,3)); % for legend
for iAtt = 1:nAtt
    m = squeeze(mean(vhDiffSTAve(:,iAtt,:),3));
%     se = squeeze(std(vhDiffSTAve(:,iAtt,:),0,3)/sqrt(nSubjects));
    se = squeeze(std(vhDiffSTAveAttDiff,0,2)/sqrt(nSubjects));
    shadedErrorBar(tt, m, se, ...
        {'Color',allColors.ampsAtt(iAtt,:)}, 1);
end
ylim(ylims)
vline(0,'k')
plot(tt([1 end]),[0 0],'k','LineWidth',1)
title('T1 & T2')
xlabel('Time (ms)')
ylabel('ITPC V-H')
legend(p1, attNames)

%% V-H, att-unatt
for iT = 1:nT
    [h p] = ttest(squeeze(vhDiffSAttDiff(:,iT,:))');
    figure
    hold on
    plot(tt, squeeze(vhDiffSAttDiff(:,iT,:)),'Color',allColors.ampsPA(1,:))
    plot(tt, squeeze(mean(vhDiffSAttDiff(:,iT,:),3)),'Color',allColors.ampsPA(1,:),'LineWidth',4)
    plot(tt, h,'Color',[.3 .3 .3])
    plot(twin,[0 0],'k')
    ylim([-.15 .15])
    vline(0,'k')
    title(sprintf('%s, V-H att - V-H unatt',targetNames{iT}))
end

%% V-H, att and unatt overlaid, T1 & T2 average
figure
hold on
for iAtt = 1:nAtt
    plot(tt, squeeze(vhDiffSTAve(:,iAtt,:)), 'Color',allColors.ampsAtt(iAtt,:))
end
for iAtt = 1:nAtt
    p1(iAtt) = plot(tt, squeeze(mean(vhDiffSTAve(:,iAtt,:),3)), 'Color',allColors.ampsAtt(iAtt,:),'LineWidth',4);
end
plot(tt([1 end]),[0 0],'k')
vline(0,'k')
legend(p1, 'att','unatt')
title('T1 & T2 average, V-H')

%% V-H, att-unatt, T1 & T2 average
[h p] = ttest(vhDiffSTAveAttDiff');

figure
hold on
plot(tt, vhDiffSTAveAttDiff, 'Color',allColors.ampsPA(1,:))
plot(tt, mean(vhDiffSTAveAttDiff,2), 'Color',allColors.ampsPA(1,:),'LineWidth',4)
plot(tt, h, 'Color',[.3 .3 .3])
plot(tt([1 end]),[0 0],'k')
ylim([-.15 .15])
vline(0,'k')
title('T1 & T2 average, V-H att - V-H unatt')

%% save analysis
vhFileName = sprintf('%s/gN%d_VH_workspace.mat', analysisDir, nSubjects);

if saveAnalysis
    save(vhFileName)
end


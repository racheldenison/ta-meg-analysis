% rd_TA2_plotDecodeGroup.m

%% setup
expt = 'TA2';

analysisName = 'classAccT1T2ByCond_sp10_nt3';

switch expt
    case 'TANoise'
        exptDir = '/Local/Users/denison/Data/TANoise/MEG';
        
        subjects = {'R0817_20171212','R0817_20171213',...
            'R1187_20180105','R1187_20180108',...
            'R0983_20180111','R0983_20180112',...
            'R0898_20180112','R0898_20180116',...
            'R1021_20180208','R1021_20180212',...
            'R1103_20180213','R1103_20180215',...
            'R0959_20180219','R0959_20180306'}; % N=7 x 2 sessions TANoise
        
    case 'TA2'
        exptDir = '/Local/Users/denison/Data/TA2/MEG';
        
        subjects = {'R0817_20181120','R0890_20181121','R0959_20181128',...
            'R1103_20181121','R1187_20181119','R1373_20181128','R1452_20181119'};
end
nSubjects = numel(subjects);

saveFileName = sprintf('%s/Group/mat/gN%d_%s.mat', exptDir, nSubjects, analysisName);

saveData = 0;

%% get data
groupData = [];
for iSubject = 1:nSubjects
    sessionDir = subjects{iSubject};
    analysisFileName = sprintf('%s/%s/mat/%s.mat', exptDir, sessionDir, analysisName);
    
    load(analysisFileName)
    
    groupData(:,:,:,iSubject) = A.classAcc;
end

cueNames = A.cueNames;
targetNames = A.targetNames;
targetWindows = A.targetWindows;
sp = A.decodingOps.binSize;
nTarget = numel(targetWindows);

%% summarize
groupMean = mean(groupData,4);
groupSte = std(groupData,0,4)./sqrt(nSubjects);

groupDiff = squeeze(groupData(:,1,:,:) - groupData(:,2,:,:)); % precueT1 - precueT2
groupDiffMean = mean(groupDiff,3);
groupDiffSte = std(groupDiff,0,3)./sqrt(nSubjects);

groupDiffT1T2Ave = squeeze((groupDiff(:,1,:) - groupDiff(:,2,:))./2);
groupDiffT1T2AveMean = mean(groupDiffT1T2Ave,2);
groupDiffT1T2AveSte = std(groupDiffT1T2Ave,0,2)./sqrt(nSubjects);

%% save group data
if saveData
    save(saveFileName,'A','groupData')
end

%% colors
figure
colors = get(gca,'colororder');
close(gcf)

%% plot
ylims = [40 60];

figure
for iT = 1:nTarget
    twin = targetWindows{iT};
    times = twin(1):sp:twin(2);
    xlims = twin;
    
    subplot(2,1,iT)
    hold on
    plot(times, groupMean(:,1:2,iT),'LineWidth',1)
    for iCond=1:2
        shadedErrorBar(times, groupMean(:,iCond,iT), groupSte(:,iCond,iT), {'color', colors(iCond,:), 'LineWidth', 3}, 1)
    end
    plot(xlims,[50 50],'k')
    xlim(xlims)
    ylim(ylims)
    xlabel('time (ms)')
    title(sprintf('T%d',iT))
    if iT==1
        legend(cueNames(1:2))
    else
        ylabel('classification accuracy (%)')
    end
end


times = 0:sp:targetWindows{1}(2)-targetWindows{1}(1);
xlims = times([1 end]);
ylims = [-8 8];

figure
hold on
plot(times, groupDiffMean,'LineWidth',1)
for iT=1:nTarget
    shadedErrorBar(times, groupDiffMean(:,iT), groupDiffSte(:,iT), {'color', colors(iT,:), 'LineWidth', 3}, 1)
end
plot(xlims,[0 0],'k')
xlim(xlims)
ylim(ylims)
xlabel('time (ms)')
ylabel({'classification accuracy (%)','precue T1 - precue T2'})
legend('T1','T2')
 

figure
for iT = 1:nTarget
    if iT==1
        vals = groupDiffMean(:,iT);
    else
        vals = -groupDiffMean(:,iT);
    end
    subplot(2,1,iT)
    hold on
    shadedErrorBar(times, vals, groupDiffSte(:,iT), {'color', [0 0 0], 'LineWidth', 3}, 1)
    plot(xlims,[0 0],'k')
    xlim(xlims)
    ylim(ylims)
    title(targetNames{iT})
end
xlabel('time (ms)')
ylabel({'classification accuracy (%)','cued - uncued'})


figure
hold on
shadedErrorBar(times, groupDiffT1T2AveMean, groupDiffT1T2AveSte, {'color', [0 0 0], 'LineWidth', 3}, 1)
plot(xlims,[0 0],'k')
xlim(xlims)
xlabel('time (ms)')
ylabel({'classification accuracy (%)','cued - uncued'})
title('T1 & T2')

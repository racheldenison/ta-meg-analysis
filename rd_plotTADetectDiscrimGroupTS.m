function rd_plotTADetectDiscrimGroupTS(A, measure, groupData, groupMean, groupSte)

%% setup
plotOrder = [1 5 3 7 2 6 4 8 9];

tf9FigPos = [0 250 1280 580];
tf3FigPos = [200 475 1000 275];

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
        ylims = [-.5 2.5];
        diffYlims = [-.8 .8];
    case 'h'
        ylims = [-5 40];
        diffYlims = [-5 5];
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

%% calculate attT2-T1
valsDiff = diff(groupData.ampsAtt);
valsDiffMean = squeeze(mean(valsDiff,3));
valsDiffSte = squeeze(std(valsDiff,0,3)./sqrt(nSubjects));

%% indiv subjects ts
for iF = 1:nFields
    fieldName = fieldNames{iF};
    
    vals = groupData.(fieldName);
    valsMean = groupMean.(fieldName);
    colors = allColors.(fieldName);
    
    figure
    set(gcf,'Position',tf9FigPos)
    for iSubject = 1:nSubjects
        subplot(nrows,ncols,iSubject)
        hold on
        if strcmp(fieldName, 'amps') % condition dimension varies (ick)
            for iCond=1:size(valsMean,2)
                plot(t, vals(:,iCond,iSubject), 'color', colors(iCond,:))
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
    end
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
    ylim(diffYlims)
    for iEv = 1:numel(eventTimes)
        vline(eventTimes(iEv),'color','k','LineStyle',':');
    end
    if iSubject==1
        xlabel('time (ms)')
        ylabel('amplitude difference (T2-T1)')
    end
end

%% group
for iF = 1:nFields
    fieldName = fieldNames{iF};
    
    valsMean = groupMean.(fieldName);
    valsSte = groupSte.(fieldName);
    colors = allColors.(fieldName);
    
    if strcmp(fieldName, 'amps')
        valsMean = valsMean';
        valsSte = valsSte';
    end
    
    figure
    hold on
    for iCond=1:size(valsMean,1)
        shadedErrorBar(t, valsMean(iCond,:), valsSte(iCond,:), {'color', colors(iCond,:), 'LineWidth', 3}, 1)
    end
    plot(xlims, [0 0], 'k')
    for iEv = 1:numel(eventTimes)
        vline(eventTimes(iEv),'color','k','LineStyle',':');
    end
    xlabel('time (ms)')
    ylabel('amplitude')
end

%% group attT2-attT1 with ste error bars
valsMean = groupMean.ampsAtt;
colors = allColors.ampsAtt;
figure
hold on
for iCond=1:2
    shadedErrorBar(t, valsMean(iCond,:), valsDiffSte, {'color', colors(iCond,:), 'LineWidth', 3}, 1)
end
plot(xlims, [0 0], 'k')
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'color','k','LineStyle',':');
end
xlabel('time (ms)')
ylabel('amplitude')

figure
hold on
shadedErrorBar(t, valsDiffMean, valsDiffSte, 'k', 1)
plot(xlims, [0 0], 'k')
ylim(diffYlims)
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'color','k','LineStyle',':');
end
xlabel('time (ms)')
ylabel('amplitude difference (T2-T1)')


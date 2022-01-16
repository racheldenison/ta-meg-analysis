% rd_plotPeaks2.m

exptDir = pathToTANoise('MEG');
dataDir = sprintf('%s/Group/mat', exptDir);
% dataDir = '/Local/Users/denison/Data/TANoise/MEG/Group/mat';
% dataDir = '~/Downloads';

collapseSessions = 0;

load(sprintf('%s/gN10_peaks_20Hz.mat', dataDir))

if collapseSessions
    load(sprintf('%s/gN10_itpcAtt_20Hz.mat', dataDir))
else
    load(sprintf('%s/gN10_itpcAtt_20Hz_bySession.mat', dataDir))
end

t = t(1:6701);

tsWin = [1100 1750];
tsidx = find(t==tsWin(1)):find(t==tsWin(2));

peakWin = 100;

evTimes = [0 1050 1350 2300];
tslims = [-100 2400];

nSubjects = numel(subjects);

subjectFactors = [1 -1 -1 1 0 1 -1 1 -1 -1]; % 1 = pos, -1 = neg, 0 = none
% subjectFactors = [1 0 0 1 1 1 1];
% subjectFactors = [1 0 0 1 0 1 1 1 1 1]; % 15 Hz
% subjectFactors = [-1 1 0 1 0 1 1 0 1 1]; % 40 Hz

exampleSubjects = [1 7];

%% define data
data = itpc;

dataMean = mean(data,3);

%% normalize data (currently only used to calculate baseline for other analyses)
% baseline
btwin = [500 1000];
btidx = find(t==btwin(1)):find(t==btwin(2));
baseline = squeeze(mean(mean(data(btidx,:,:,:),1),2));

%% plot the data
figure
hold on
plot(t, dataMean)
% xlim(tslims)
for iEv = 1:numel(evTimes)
    vline(evTimes(iEv),'k')
end
xlabel('Time (ms)')
ylabel('ITPC')
legend('Precue T1','Precue T2')
legend boxoff

figure
hold on
colors = get(gca,'ColorOrder');
colors = repmat(colors,2,1);
plot(t, squeeze(data(:,1,:)),'color',colors(1,:)) % cue T1
plot(t, squeeze(data(:,2,:)),'color',colors(2,:)) % cue T2
xlabel('Time (ms)')
ylabel('ITPC')

figure
hold on
plot(t, squeeze(data(:,1,:)),'LineStyle','-') % cue T1
ax = gca;
ax.ColorOrderIndex = 1;
plot(t, squeeze(data(:,2,:)),'LineStyle','--') % cue T2
xlim([1000 2000])
xlabel('Time (ms)')
ylabel('ITPC')

figure('Position',[100 20 350 600])
for iSubject = 1:nSubjects
    subplot(nSubjects,1,iSubject)
    hold on
    plot(t, squeeze(data(:,1,iSubject)),'color',colors(1,:)) % cue T1
    plot(t, squeeze(data(:,2,iSubject)),'color',colors(2,:)) % cue T2
%     xlim([500 2000])
    ylim([0 .7])
    if iSubject==nSubjects
        xlabel('Time (ms)')
        ylabel('ITPC')
    else
        set(gca,'XTick',[])
    end
end
ax = get(gcf,'Children');
ax(7).YLim = [.3 .8];
ax(6).YLim = [.3 .8];

%% plot example subjects
ylims = [.3 .7; .1 .5];

figure
for iS = 1:numel(exampleSubjects)
    subplot(numel(exampleSubjects),1,iS)
    hold on
    plot(t, data(:,:,exampleSubjects(iS)))
    xlim(tslims)
    ylim(ylims(iS,:))
    for iEv = 1:numel(evTimes)
        vline(evTimes(iEv),'k')
    end
end
xlabel('Time (ms)')
ylabel('ITPC')
legend('precue T1','precue T2')
legend boxoff


%% plot the flipped data
for iS = 1:nSubjects
    flipData(:,:,iS) = data(:,:,iS)*subjectFactors(iS);
end
flipData(:,:,subjectFactors==0) = NaN;

flipDataMean = nanmean(flipData,3);

figure
hold on
colors = get(gca,'ColorOrder');
colors = repmat(colors,2,1);
plot(t, squeeze(flipData(:,1,:)),'color',colors(1,:)) % cue T1
plot(t, squeeze(flipData(:,2,:)),'color',colors(2,:)) % cue T2
xlabel('Time (ms)')
ylabel('ITPC flipped')

figure
hold on
plot(t, flipDataMean)
for iEv = 1:numel(evTimes)
    vline(evTimes(iEv),'k')
end
xlabel('Time (ms)')
ylabel('ITPC flipped')
xlim(tslims)
legend('precue T1','precue T2')
legend boxoff

%% flipped with baseline
% baseline
btwin = [500 1000];
btidx = find(t==btwin(1)):find(t==btwin(2));
b = mean(mean(flipData(btidx,:,:),1),2);
flipDataB = flipData - repmat(b,size(data,1),size(data,2),1);

flipDataBMean = nanmean(flipDataB,3);
flipDataBSte = nanstd(flipDataB,0,3)./sqrt(nnz(subjectFactors~=0));

figure
hold on
colors = get(gca,'ColorOrder');
colors = repmat(colors,2,1);
plot(t, squeeze(flipDataB(:,1,:)),'color',colors(1,:)) % cue T1
plot(t, squeeze(flipDataB(:,2,:)),'color',colors(2,:)) % cue T2
xlabel('Time (ms)')
ylabel('ITPC flipped and baseline-corrected')

figure
hold on
plot(t, flipDataBMean)
% shadedErrorBar(t, flipDataBMean(:,1), flipDataBSte(:,1),{'color',colors(1,:)},1)
% shadedErrorBar(t, flipDataBMean(:,2), flipDataBSte(:,2),{'color',colors(2,:)},1)
for iEv = 1:numel(evTimes)
    vline(evTimes(iEv),'k')
end
xlabel('Time (ms)')
ylabel('ITPC flipped and baseline-corrected')
xlim(tslims)
% ylim([-.1 .2])
legend('precue T1','precue T2')
legend boxoff

%% find the peaks in the right timeseries window
measures = {'peaksPos','peaksNeg'};
nM = numel(measures);
for iM = 1:nM
    m = measures{iM};
    name = sprintf('%sInWin',m);
    for iSubject = 1:nSubjects
        p = peaks(iSubject).(m);
        p(p<tsWin(1) | p>tsWin(2)) = [];
        peaks(iSubject).(name) = p;
    end
end

%% plot the peak times
measures = {'peaksPosInWin','peaksNegInWin'};
nM = numel(measures);
figure
for iM = 1:nM
    m = measures{iM};
    subplot(nM,1,iM)
    hold on
    for iSubject = 1:nSubjects
        if ~isempty(peaks(iSubject).(m))
            y = repmat(iSubject,1,length(peaks(iSubject).(m)));
            plot(peaks(iSubject).(m),y,'Color',colors(iSubject,:))
            plot(peaks(iSubject).(m),y,'.','MarkerSize',30,'Color',colors(iSubject,:))
        end
    end
    xlim(tsWin)
    xlabel('time (ms)')
    ylabel('observer')
    title(m)
end

%% get the final peaks
for iSubject = 1:nSubjects
    f = subjectFactors(iSubject);
    switch f
        case 1
            peaks(iSubject).peaksSelected = peaks(iSubject).peaksPosInWin(1:2);
        case -1
            peaks(iSubject).peaksSelected = peaks(iSubject).peaksNegInWin(1:2);
        case 0
            peaks(iSubject).peaksSelected = [];
        otherwise
            error('subject factor not found')
    end
end

%% plot final peaks
figure
hold on
for iSubject = 1:nSubjects
    if ~isempty(peaks(iSubject).peaksSelected)
        plot(peaks(iSubject).peaksSelected,[iSubject iSubject],'Color',colors(iSubject,:))
        plot(peaks(iSubject).peaksSelected,[iSubject iSubject],'.','MarkerSize',30,'Color',colors(iSubject,:))
    end
end
% xlim(tsWin)
xlim([1000 1600])
ylim([0 nSubjects+1])
xlabel('Time (ms)')
ylabel('Observer')
title('Selected peaks')
set(gca,'YTick',1:nSubjects)

% add target time lines, adjusted by 50 ms to account for pd delay
vline(evTimes(2),'--k')
vline(evTimes(3),'--k')


%% extract timeseries in the window around the peak
tsdata = flipDataB;

peakData = [];
for iSubject = 1:nSubjects
%     f = subjectFactors(iSubject);
    f = 1;
    p = peaks(iSubject).peaksSelected;
    for iP = 1:numel(p)
        tidx = find(t==p(iP)-peakWin/2):find(t==p(iP)+peakWin/2);
        if collapseSessions
            vals = tsdata(tidx,:,iSubject);
            peakData(:,:,iP,iSubject) = vals*f; % [time cond peak subject]
        else
            vals = tsdata(tidx,:,iSubject,:);
            peakData(:,:,iP,iSubject,:) = vals*f; % [time cond peak subject session]
        end
    end
end
peakData(:,:,:,subjectFactors==0,:) = NaN;

%% plot the peak timeseries
if collapseSessions
    figure
    for iP = 1:2
        subplot(1,2,iP)
        hold on
        for iSubject = 1:nSubjects
            plot(peakData(:,1,iP,iSubject),'color',colors(1,:))
            plot(peakData(:,2,iP,iSubject),'color',colors(2,:))
        end
        xlim([1 size(peakData,1)])
    end
else
    figure
    for iSession = 1:2
        for iP = 1:2
            subplot(2,2,iP+(iSession-1)*2)
            hold on
            for iSubject = 1:nSubjects
                plot(peakData(:,1,iP,iSubject,iSession),'color',colors(1,:))
                plot(peakData(:,2,iP,iSubject,iSession),'color',colors(2,:))
            end
            xlim([1 size(peakData,1)])
            if iSession==1
                title(sprintf('Peak %d', iP))
            end
            if iP==1
                ylabel(sprintf('Session %d', iSession))
            end
        end
    end
end

%% summarize
peakDataAve = squeeze(mean(peakData,1));
peakMean = nanmean(peakDataAve,3);
peakSte = nanstd(peakDataAve,0,3)./sqrt(nnz(subjectFactors~=0));
if collapseSessions
    peakDataDiff = squeeze(peakDataAve(1,:,:) - peakDataAve(2,:,:));
    peakDiffSte = nanstd(peakDataDiff,0,2)./sqrt(nnz(subjectFactors~=0));
else
%     save(sprintf('%s/gN10_itpcAttPeakVals_20Hz_bySession.mat', dataDir), 'peakData','peakDataAve')
end

%% plots
if collapseSessions
    figure
    bar(peakMean')
    xlabel('peak')
    ylabel('ITPC (avearage after flipping)')
    legend('precue T1','precue T2')
    
    figure('Position',[150 200 350 450])
    hold on
    for iCue = 1:2
        errorbar(1:2, peakMean(iCue,:), peakDiffSte, '.', 'MarkerSize', 30)
%         errorbar(1:2, peakMean(iCue,:), peakSte(iCue,:), '.', 'MarkerSize', 30)
    end
    ax = gca;
    jitterx(ax);
    ax.XTick = [1 2];
    ylim([.05 .1])
    ax.YTick = .05:.01:1;
    xlabel('Peak')
    ylabel('ITPC (average after flipping and baseline correction)')
    legend('precue T1','precue T2','Location','best')
    legend boxoff
    
    figure
    for iP = 1:2
        subplot(1,2,iP)
        hold on
        plot(squeeze(peakDataAve(1,iP,:)),squeeze(peakDataAve(2,iP,:)),'.','MarkerSize', 30)
        plot([-1 1],[-1 1],'k')
        axis square
        xlabel('precue T1')
        ylabel('precue T2')
    end
    
    figure
    bar(peakDataDiff')
    xlabel('Observer')
    ylabel('ITPC, precue T1 - precue T2')
    legend('peak 1','peak 2')
end

%% plot baseline for uppers and downers
bup = baseline(subjectFactors==1,:);
bdown = baseline(subjectFactors==-1,:);

figure
hold on
plot(ones(size(bup(:))), bup(:), 'o', 'MarkerSize',12, 'MarkerFaceColor','w','MarkerEdgeColor','k')
plot(ones(size(bdown(:)))*2, bdown(:), 'o', 'MarkerSize',16, 'MarkerFaceColor','k','MarkerEdgeColor','w')
xlim([0 3])
set(gca,'XTick',[1 2],'XTickLabel',{'Upward','Downward'})
ylabel('Baseline ITPC')
xlabel('Peak direction')

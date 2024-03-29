% rd_plotPeaks.m

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

tsWin = [1050 1650]; %[1100 1750];
tsidx = find(t==tsWin(1)):find(t==tsWin(2));

peakWin = 100;

nSubjects = numel(subjects);

subjectFactors = [1 -1 -1 1 0 1 -1 1 -1 -1]; % 1 = pos, -1 = neg, 0 = none
% subjectFactors = [1 0 0 1 1 1 1];
% subjectFactors = [1 0 0 1 0 1 1 1 1 1]; % 15 Hz
% subjectFactors = [-1 1 0 1 0 1 1 0 1 1]; % 40 Hz

%% define data
data = itpc;

%% plot the data
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
xlim(tsWin)
xlabel('time (ms)')
ylabel('observer')
title('selected peaks')

%% extract timeseries in the window around the peak
peakData = [];
for iSubject = 1:nSubjects
    f = subjectFactors(iSubject);
    p = peaks(iSubject).peaksSelected;
    for iP = 1:numel(p)
        tidx = find(t==p(iP)-peakWin/2):find(t==p(iP)+peakWin/2);
        if collapseSessions
            vals = data(tidx,:,iSubject);
            peakData(:,:,iP,iSubject) = vals*f; % [time cond peak subject]
        else
            vals = data(tidx,:,iSubject,:);
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
if collapseSessions
    peakDataDiff = squeeze(peakDataAve(1,:,:) - peakDataAve(2,:,:));
    peakDataDiffSte = nanstd(peakDataDiff,0,2)./sqrt(nnz(subjectFactors~=0));
else
    save(sprintf('%s/gN10_itpcAttPeakVals_20Hz_bySession.mat', dataDir), 'peakData','peakDataAve')
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
        errorbar(1:2, peakMean(iCue,:), peakDataDiffSte, '.', 'MarkerSize', 30)
    end
    ax = gca;
    jitterx(ax);
    ax.XTick = [1 2];
    xlabel('Peak')
    ylabel('ITPC (average after flipping)')
    legend('precue T1','precue T2','Location','best')
    
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
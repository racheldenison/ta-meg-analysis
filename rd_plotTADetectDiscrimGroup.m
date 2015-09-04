% rd_plotTADetectDiscrimGroup.m

%% Setup
exptDir = '/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG';
analStr = 'ebi_ft'; % '', 'ebi', etc.
ssvefFreq = 30;
nTopChannels = 1; % 1, 5, etc.

% subjects = {'R0817_20150504', 'R0973_20150727', 'R0974_20150728', ...
%     'R0861_20150813', 'R0983_20150813', 'R0898_20150828'};
subjects = {'R0817_20150504', 'R0973_20150727', ...
    'R0861_20150813', 'R0983_20150813', 'R0898_20150828'};

nSubjects = numel(subjects);

tstart = -500; % ms 
tstop = 3600; % ms 
t = tstart:tstop;
evTimes = [0 500 1500 2100 3100];

%% Get data
for iSubject = 1:nSubjects
    subject = subjects{iSubject};
    
    sessionDir = subject;
    dataDir = sprintf('%s/%s', exptDir, sessionDir);
    matDir = sprintf('%s/mat', dataDir);
    
    analysisFile = dir(sprintf('%s/analysis_*_topChannels%d_%dHz.mat', matDir, nTopChannels, ssvefFreq));

    if numel(analysisFile)==1
        load(sprintf('%s/%s', matDir, analysisFile.name))
    else
        error('too many or too few matching analysis files')
    end
    
    wAmps(:,:,iSubject) = A.wAmpsAtt;
end

%% Calculate
wAmpsDiff = diff(wAmps);

%% summary stats
wAmpsMean = squeeze(mean(wAmps,3));
wAmpsSte = squeeze(std(wAmps,0,3)./sqrt(nSubjects));

wAmpsDiffMean = squeeze(mean(wAmpsDiff,3));
wAmpsDiffSte = squeeze(std(wAmpsDiff,0,3)./sqrt(nSubjects));

%% plot figs
colors = {'b','r'};
xlims = [t(1) t(end)];
ylims = [-.5 2.5];
diffLims = [-1 1];

% indiv subjects
figure
set(gcf,'Position',[30 450 1200 250])
for iSubject = 1:nSubjects
    subplot(1,nSubjects,iSubject)
    hold on
    for i=1:2
        plot(t, wAmps(i,:,iSubject), colors{i})
    end
    plot(xlims, [0 0], 'k')
    xlim(xlims)
    ylim(ylims)
    for iEv = 1:numel(evTimes)
        vline(evTimes(iEv),'color','k','LineStyle',':');
    end
end

figure
set(gcf,'Position',[30 450 1200 250])
for iSubject = 1:nSubjects
    subplot(1,nSubjects,iSubject)
    hold on
    plot(t, wAmpsDiff(1,:,iSubject), 'k', 'LineWidth', 2)
    plot(xlims, [0 0], 'k')
    xlim(xlims)
    ylim(diffLims)
    for iEv = 1:numel(evTimes)
        vline(evTimes(iEv),'color','k','LineStyle',':');
    end
end

% group
figure
hold on
for i=1:2
    shadedErrorBar(t, wAmpsMean(i,:), wAmpsSte(i,:), colors{i}, 1)
end
plot(xlims, [0 0], 'k')
for iEv = 1:numel(evTimes)
    vline(evTimes(iEv),'color','k','LineStyle',':');
end

figure
hold on
for i=1:2
    shadedErrorBar(t, wAmpsMean(i,:), wAmpsDiffSte, colors{i}, 1)
end
plot(xlims, [0 0], 'k')
for iEv = 1:numel(evTimes)
    vline(evTimes(iEv),'color','k','LineStyle',':');
end

figure
hold on
shadedErrorBar(t, wAmpsDiffMean, wAmpsDiffSte, 'k', 1)
plot(xlims, [0 0], 'k')
ylim([-.5 .5])
for iEv = 1:numel(evTimes)
    vline(evTimes(iEv),'color','k','LineStyle',':');
end


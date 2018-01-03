% rd_plotMultiSessionData.m

exptType = 'TANoise';
% a1 = load('/Local/Users/denison/Data/TANoise/MEG/R0817_20171212/mat/analysis_R0817_TANoise_12.12.17_ebi_ft_topChannels5_allTrials_20Hz.mat');
% a2 = load('/Local/Users/denison/Data/TANoise/MEG/R0817_20171213/mat/analysis_R0817_TANoise_12.13.17_ebi_ft_topChannels5_allTrials_20Hz.mat');

a1 = load('/Local/Users/denison/Data/TANoise/MEG/R0817_20171212/mat/analysis_singleTrials_R0817_TANoise_12.12.17_ebi_ft_topChannels5_allTrials_20Hz.mat');
a2 = load('/Local/Users/denison/Data/TANoise/MEG/R0817_20171213/mat/analysis_singleTrials_R0817_TANoise_12.13.17_ebi_ft_topChannels5_allTrials_20Hz.mat');

A(1) = a1.A;
A(2) = a2.A;
nA = numel(A);

trigNames = A(1).trigNames;
nTrigs = numel(trigNames);
t = A(1).t;
eventTimes = A(1).eventTimes;
attNames = A(1).attNames;

%% Plotting setup
plotOrder = [1 5 3 7 2 6 4 8 9];
extendedMap = flipud(lbmap(nTrigs-1+4,'RedBlue'));
selectedMap = extendedMap([1:(nTrigs-1)/2 (end-(nTrigs-1)/2)+1:end],:);
trigColors = [selectedMap; 0 0 0];
trigBlue = mean(selectedMap(1:(nTrigs-1)/2,:));
trigRed = mean(selectedMap((end-(nTrigs-1)/2)+1:end,:));

tsFigPos = [0 500 1250 375];
% ts2FigPos = [0 500 1100 600];
% ts3FigPos = [0 500 1100 900];
% condFigPos = [250 300 750 650];
% tf9FigPos = [0 250 1280 580];
tf3FigPos = [200 475 980 330];

set(0,'defaultLineLineWidth',1)

switch exptType
    case 'TADetectDiscrim'
        xtickint = 50;
    case 'TAContrast'
        xtickint = 100;
    case 'TANoise'
        xtickint = 100;
end


%% wAmps
vals = [];
for iA = 1:nA
    vals(:,:,iA) = A(iA).wAmps;
end
wAmps = nanmean(vals, 3);

fH = [];
fH(1) = figure;
set(gcf,'Position',tsFigPos)
set(gca,'ColorOrder',trigColors)
hold all
% plot(t, wAmps(:,plotOrder))
plot(t, wAmps(:,end), 'k') % blank
plot(t, nanmean(wAmps(:,plotOrder(1:(nTrigs-1)/2)),2),'color',trigBlue,'LineWidth',4)
plot(t, nanmean(wAmps(:,plotOrder(end-(nTrigs-1)/2):end-1),2),'color',trigRed,'LineWidth',4)
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
% legend(trigNames(plotOrder))
legend('blank','att T1','att T2')
xlabel('time (ms)')
ylabel('wavelet amp')
% title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels) wstrt])

%% time freq single
vals = [];
for iA = 1:nA
    vals(:,:,:,iA) = A(iA).stfAmpsAtt;
end
tfSingleAmpsAtt = nanmean(vals, 4);

% figures
toi = A(1).stfToi;
foi = A(2).stfFoi;
ytick = 10:10:numel(foi);
xtick = 51:xtickint:numel(toi);
clims = [0 70];
diffClims = [-10 10];

fH(2) = figure;
set(gcf,'Position',tf3FigPos)
attNames = {'attT1','attT2'};
for iAtt = 1:size(tfSingleAmpsAtt,3)
    subplot(1,3,iAtt)
    imagesc(tfSingleAmpsAtt(:,:,iAtt),clims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    xlabel('time (s)')
    ylabel('frequency (Hz)')
    title(attNames{iAtt})
end
subplot(1,3,3)
imagesc(tfSingleAmpsAtt(:,:,2)-tfSingleAmpsAtt(:,:,1),diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (s)')
ylabel('frequency (Hz)')
title('attT2 - attT1')
% rd_supertitle(['channel' sprintf(' %d', channels) wstrt]);
rd_raiseAxis(gca);

%% wAmps single
wAmpsAtt = cat(2, A(1).wAmpsAtt, A(2).wAmpsAtt);

fH(3) = figure;
set(gcf,'Position',tsFigPos)
hold on
plot(t, nanmean(wAmpsAtt(:,:,1),2),'color',trigBlue,'LineWidth',4)
plot(t, nanmean(wAmpsAtt(:,:,2),2),'color',trigRed,'LineWidth',4)
legend(attNames)
[~, emp, err] = rd_bootstrapCI(wAmpsAtt(:,:,1)');
shadedErrorBar(t, emp, err, {'color',trigBlue,'LineWidth',4}, 1)
[~, emp, err] = rd_bootstrapCI(wAmpsAtt(:,:,2)');
shadedErrorBar(t, emp, err, {'color',trigRed,'LineWidth',4}, 1)
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
xlabel('time (ms)')
ylabel('single trial wavelet amp')

%% itpc - NB plotting means, not recomputing ITPC!
vals = [];
for iA = 1:nA
    vals(:,:,iA) = A(iA).wITPCAtt;
end
wITPCAtt = nanmean(vals, 3);

fH(4) = figure;
set(gcf,'Position',tsFigPos)
hold on
plot(t, wITPCAtt(:,1),'color',trigBlue,'LineWidth',4)
plot(t, wITPCAtt(:,2),'color',trigRed,'LineWidth',4)
legend(attNames)
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
xlabel('time (ms)')
ylabel('wavelet itpc')


function rd_plotTADetectDiscrimGroupTS(A, groupMean, saveFigs, figDir, figStr)

%% args
if nargin<3
    saveFigs = 0;
end
if nargin<5
    figStr = '';
    if saveFigs==1
        error('If youre saving figs, you should specify a figStr')
    end
end

figTitle = und2space(figStr);

%% setup
plotOrder = [1 5 3 7 2 6 4 8 9];

ts2FigPos = [0 500 1100 600];
ts3FigPos = [0 500 1100 900];

eventTimes = A.eventTimes;
trigNames = A.trigNames;

t = A.t;
Fs = A.Fs;
tidx1 = find(t==eventTimes(2));
tidx2 = find(t==eventTimes(5))-1;
nfft = numel(tidx1:tidx2);
f = Fs/2*linspace(0,1,nfft/2+1);
targetF = A.targetF;
targetWindow = A.targetWindow;

nTrigs = numel(trigNames);

% colors
extendedMap = flipud(lbmap(nTrigs-1+4,'RedBlue'));
selectedMap = extendedMap([1:(nTrigs-1)/2 (end-(nTrigs-1)/2)+1:end],:);
trigColors = [selectedMap; 0 0 0];

set(0,'defaultLineLineWidth',1)

%% Time series and FFT
figure
set(gcf,'Position',ts2FigPos)

% mean across channels
trigMean = mean(groupMean.tsAmps,2);
amps = mean(groupMean.fAmps,2);

% mean across attention conditions
tsAttT1 = mean(trigMean(:,plotOrder(1:4)),2);
tsAttT2 = mean(trigMean(:,plotOrder(5:8)),2);

% time
subplot(3,1,1)
set(gca,'ColorOrder',trigColors)
hold all
plot(t, squeeze(trigMean(:,1,plotOrder)))
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
xlim([t(1) t(end)])
xlabel('time (ms)')
ylabel('amplitude')
title(figTitle)

% frequency
subplot(3,1,2)
set(gca,'ColorOrder',trigColors)
hold all
plot(repmat(f',1,nTrigs), squeeze(amps(:,1,plotOrder)))
xlim([1 200])
ylim([0 20])
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
legend(trigNames(plotOrder))

subplot(3,1,3)
set(gca,'ColorOrder',[.66 .5 .78; trigColors(end,:)])
hold all
plot(f, squeeze(mean(amps(:,1,1:end-1),3)))
plot(f, squeeze(amps(:,1,end)))
xlim([1 200])
ylim([0 20])
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
legend('stim average','blank')

if saveFigs
    figPrefix = sprintf('%s_plot', figStr);
    rd_saveAllFigs(gcf, {'tsFFT'}, figPrefix, figDir)
end

%% Target present/absent
names = {'target present','target absent'};
targetPA = groupMean.targetPA;
targetPADiff = groupMean.targetPADiff;
targetAmps = groupMean.targetPADiffAmps;

fH = [];
fH(1) = figure;
set(gcf,'Position',ts2FigPos)
for iPA = 1:2
    subplot(2,1,iPA)
    plot(targetWindow(1):targetWindow(2), targetPA(:,:,iPA))
    vline(0,'k');
    if iPA==1
        legend('xx1','xx2','xo1','ox2')
    end
    xlabel('time (ms)')
    ylabel('amplitude')
    title(names{iPA})
end
rd_supertitle(figTitle)
rd_raiseAxis(gca);

fH(2) = figure;
set(gcf,'Position',ts3FigPos)
subplot(3,1,1)
plot(targetWindow(1):targetWindow(2), targetPADiff)
xlabel('time (ms)')
ylabel('\Delta amplitude')
title('target present - absent')
subplot(3,1,2)
plot(targetWindow(1):targetWindow(2), mean(targetPADiff,1), 'k');
xlabel('time (ms)')
ylabel('\Delta amplitude')
subplot(3,1,3)
plot(targetF, targetAmps)
xlim([0 150])
xlabel('frequency (Hz)')
ylabel('\Delta amplitude')
rd_supertitle(figTitle)
rd_raiseAxis(gca);

fH(3) = figure;
hold on
plot(t, tsAttT1, 'color', mean(trigColors(1:4,:)))
plot(t, tsAttT2, 'color', mean(trigColors(5:8,:)))
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
xlim([t(1) t(end)])
xlabel('time (ms)')
ylabel('amplitude')
title(figTitle)
legend('attT1','attT2')


if saveFigs
    figPrefix = sprintf('%s_plot', figStr);
    rd_saveAllFigs(fH, {'targetPATrialAve','targetPATrialAveDiff'}, figPrefix, figDir)
end

% rd_TADetectDiscrimSSVEF2.m

%% Setup
exptDir = '/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG';
sessionDir = 'R0974_20150728';
fileBase = 'R0974_TADeDi_7.28.15';
analStr = 'ebi'; % '', 'ebi', etc.
ssvefFreq = 40;
topChannels = 1:5; % 1, 1:5, etc.

dataDir = sprintf('%s/%s', exptDir, sessionDir);
matDir = sprintf('%s/mat', dataDir);

switch analStr
    case ''
        savename = sprintf('%s/%s_ssvef_workspace.mat', matDir, fileBase);
        channelsFileName = sprintf('%s/channels_%dHz.mat', matDir, ssvefFreq);
        analysisFileName = sprintf('%s/analysis_%s_topChannels%d_%dHz.mat', matDir, fileBase, numel(topChannels), ssvefFreq);
    otherwise
        savename = sprintf('%s/%s_%s_ssvef_workspace.mat', matDir, fileBase, analStr);
        channelsFileName = sprintf('%s/channels_%dHz_%s.mat', matDir, ssvefFreq, analStr);
        analysisFileName = sprintf('%s/analysis_%s_%s_topChannels%d_%dHz.mat', matDir, fileBase, analStr, numel(topChannels), ssvefFreq);
end

% eventTimes = [0 500 1500 2100 3100];

%% Get the data
load(savename)

%% Settings after loading the data
saveAnalysis = 1;
saveFigs = 1;

excludeTrialsFt = 1;
excludeSaturatedEpochs = 0;

load(channelsFileName,'channelsRanked');
channels = channelsRanked(topChannels);

%% Store settings for this analysis
A.fileBase = fileBase;
A.analStr = analStr;
A.excludeTrialsFt = excludeTrialsFt;
A.excludeSaturatedEpochs = excludeSaturatedEpochs;
A.ssvefFreq = ssvefFreq;
A.channels = channels;
A.Fs = Fs;
A.t = t;
A.eventTimes = eventTimes;
A.trigNames = trigNames;

%% Baseline
% baselinePeriod = -500:0;
baselinePeriod = t;
inBaseline = ismember(t,baselinePeriod);
baselineDC = mean(trigData(inBaseline,:,:),1);
baselineTSeries = repmat(baselineDC,[size(trigData,1),1,1]);

% trigData0 = trigData;
trigData = trigData-baselineTSeries;

%% Excluded saturated channel epochs
if excludeSaturatedEpochs
    load([matDir '/saturated_channel_epochs.mat'])
    trigData(:,saturatedChannelEpochs) = NaN;
end

%% Exclude trials manually rejected with ft
if excludeTrialsFt
    % load trials_rejected variable from ft manual rejection
    load([matDir '/trials_rejected.mat'])
    trigData(:,:,trials_rejected) = NaN;
    
    % update figDir
    figDir = [figDir '_ft'];
    
    % update analysis file
    switch analStr
        case ''
            analysisFileName = sprintf('%s/analysis_%s_ft_topChannels%d_%dHz.mat', matDir, fileBase, numel(topChannels), ssvefFreq);
        otherwise
            analysisFileName = sprintf('%s/analysis_%s_%s_ft_topChannels%d_%dHz.mat', matDir, fileBase, analStr, numel(topChannels), ssvefFreq);
    end
end

%% Make figDir if needed
if ~exist(figDir,'dir') && saveFigs
    mkdir(figDir)
end

%% Organize trials into conditions
cueCondIdx = strcmp(behav.responseData_labels, 'cue condition');
t1CondIdx = strcmp(behav.responseData_labels, 'target type T1');
t2CondIdx = strcmp(behav.responseData_labels, 'target type T2');

blankCond = 1;
cueConds = {[2 3], [4 5]}; % cue T1, cue T2
t1Conds = {[1 2], 0}; % present, absent
t2Conds = {[1 2], 0}; % present, absent

condData = [];
for iCue = 1:numel(cueConds)
    vals = cueConds{iCue};
    wCue = [];
    for iEl = 1:numel(vals)
        wCue(:,iEl) = behav.responseData_all(:,cueCondIdx) == vals(iEl);
    end
    for iT1 = 1:numel(t1Conds)
        vals = t1Conds{iT1};
        wT1 = [];
        for iEl = 1:numel(vals)
            wT1(:,iEl) = behav.responseData_all(:,t1CondIdx) == vals(iEl);
        end
        for iT2 = 1:numel(t2Conds)
            vals = t2Conds{iT2};
            wT2 = [];
            for iEl = 1:numel(vals)
                wT2(:,iEl) = behav.responseData_all(:,t2CondIdx) == vals(iEl);
            end
            
            w = sum([wCue wT1 wT2],2)==3;
            condData(:,:,:,iCue,iT1,iT2) = trigData(:,:,w);
        end
    end
end

wBlank = behav.responseData_all(:,cueCondIdx) == blankCond;
blankData = trigData(:,:,wBlank);

% mean across trials
condDataMean = squeeze(nanmean(condData,3));
blankDataMean = squeeze(nanmean(blankData,3));

% let trigMean have the conditions 1-9 in the third dimension
%%% note that here we are selecting channels!
trigMean = condDataMean(:,channels,:);
trigMean(:,:,end+1) = blankDataMean(:,channels);
nTrigs = size(trigMean,3);

A.trigMean = trigMean;

%% FFT on mean time series for each trigger type
% do the fft for each channel
nfft = 2^nextpow2(nSamples); % Next power of 2 from length of y
Y = fft(trigMean,nfft)/nSamples; % Scale by number of samples
f = Fs/2*linspace(0,1,nfft/2+1); % Fs/2 is the maximum frequency that can be measured
amps = 2*abs(Y(1:nfft/2+1,:,:)); % Multiply by 2 since only half the energy is in the positive half of the spectrum?

A.amps = amps;

%% Plotting setup
plotOrder = [1 5 3 7 2 6 4 8 9];
extendedMap = flipud(lbmap(nTrigs-1+4,'RedBlue'));
selectedMap = extendedMap([1:(nTrigs-1)/2 (end-(nTrigs-1)/2)+1:end],:);
trigColors = [selectedMap; 0 0 0];
trigBlue = mean(selectedMap(1:(nTrigs-1)/2,:));
trigRed = mean(selectedMap((end-(nTrigs-1)/2)+1:end,:));

% trigColorsPA4 = [107 76 154; 62 150 81; 57 106 177; 218 124 48]./255;
trigColorsPA4 = [.52 .37 .75; .31 .74 .40; .27 .51 .84; 1.0 .57 .22];

tsFigPos = [0 500 1250 375];
ts2FigPos = [0 500 1100 600];
ts3FigPos = [0 500 1100 900];
condFigPos = [250 300 750 650];
tf9FigPos = [0 250 1280 580];
tf3FigPos = [200 475 980 330];

set(0,'defaultLineLineWidth',1)

%% Time series and FFT
figure
set(gcf,'Position',ts2FigPos)

% time
subplot(3,1,1)
set(gca,'ColorOrder',trigColors)
hold all
plot(t, squeeze(mean(trigMean(:,:,plotOrder),2)))
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
xlim([t(1) t(end)])
xlabel('time (ms)')
ylabel('amplitude')
title(['channel' sprintf(' %d', channels)])

% frequency
subplot(3,1,2)
set(gca,'ColorOrder',trigColors)
hold all
plot(repmat(f',1,nTrigs), squeeze(mean(amps(:,:,plotOrder),2)))
xlim([1 200])
ylim([0 20])
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
legend(trigNames(plotOrder))

subplot(3,1,3)
set(gca,'ColorOrder',[.66 .5 .78; trigColors(end,:)])
hold all
plot(f, squeeze(mean(mean(amps(:,:,1:end-1),3),2)))
plot(f, squeeze(mean(amps(:,:,end),2)))
xlim([1 200])
ylim([0 20])
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
legend('stim average','blank')

if saveFigs
    if numel(channels)==1
        figPrefix = sprintf('plot_ch%d', channels);
    else
        figPrefix = ['plot_ch' sprintf('%d_', channels(1:end-1)) sprintf('%d', channels(end))];
    end
    rd_saveAllFigs(gcf, {'tsFFT'}, figPrefix, figDir)
end

%% Trial average for target present vs. absent, for a single channel
pp = mean(mean(trigMean(:,:,1:2),3),2);
pa = mean(mean(trigMean(:,:,5:6),3),2);
ap = mean(mean(trigMean(:,:,3:4),3),2);
aa = mean(mean(trigMean(:,:,7:8),3),2);

targetWindow = [-100 500];
t1Window = t>=eventTimes(3) + targetWindow(1) & t<=eventTimes(3) + targetWindow(2);
t2Window = t>=eventTimes(4) + targetWindow(1) & t<=eventTimes(4) + targetWindow(2);
targetPA(1,:,1) = pp(t1Window);
targetPA(2,:,1) = pp(t2Window);
targetPA(3,:,1) = pa(t1Window);
targetPA(4,:,1) = ap(t2Window);
targetPA(1,:,2) = aa(t1Window);
targetPA(2,:,2) = aa(t2Window);
targetPA(3,:,2) = ap(t1Window);
targetPA(4,:,2) = pa(t2Window);

targetPADiff = targetPA(:,:,1)-targetPA(:,:,2);

targetNfft = 2^nextpow2(diff(targetWindow)+1);
targetY = fft(mean(targetPADiff),targetNfft)/(diff(targetWindow)+1);
targetF = Fs/2*linspace(0,1,targetNfft/2+1);
targetAmps = 2*abs(targetY(1:targetNfft/2+1));

% store results
A.targetWindow = targetWindow;
A.targetPA = targetPA;
A.targetPADiff = targetPADiff;
A.targetF = targetF;
A.targetPADiffAmps = targetAmps;

names = {'target present','target absent'};
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
rd_supertitle(['channel' sprintf(' %d', channels)])
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
rd_supertitle(['channel' sprintf(' %d', channels)])
rd_raiseAxis(gca);

if saveFigs
    if numel(channels)==1
        figPrefix = sprintf('plot_ch%d', channels);
    else
        figPrefix = ['plot_ch' sprintf('%d_', channels(1:end-1)) sprintf('%d', channels(end))];
    end
    rd_saveAllFigs(fH, {'targetPATrialAve','targetPATrialAveDiff'}, figPrefix, figDir)
end

%% Wavelet on average across trials
switch ssvefFreq
    case 30
        width = 12; % 12 for 30 Hz, 16 for 40 Hz gives 127 ms duration, 5 Hz bandwidth
    case 40
        width = 16;
    otherwise
        error('ssvefFreq not recognized')
end
wBaselineWindow = [-500 0]; % [-300 -200];
wBaselineWindowIdx = find(t==wBaselineWindow(1)):find(t==wBaselineWindow(2));

% only frequency of interest
wAmps0 = [];
foi = ssvefFreq;
for iTrig = 1:nTrigs
    data = trigMean(:,:,iTrig)'; % channels by samples
    [spectrum,freqoi,timeoi] = ft_specest_wavelet(data, t/1000, 'freqoi', foi, 'width', width);
    specAmp = abs(squeeze(spectrum));

    if all(size(specAmp)>1) % if two-dimensional
        wAmp = specAmp;
    else
        wAmp = specAmp';
    end
    wAmpNorm = wAmp./nanmean(nanmean(wAmp(:,wBaselineWindowIdx)))-1;
    wAmps0(:,:,iTrig) = wAmpNorm';
end
wAmps = squeeze(mean(wAmps0,2)); % mean across channels

% attT1T2 means
wAmpsAtt(1,:) = mean(wAmps(:,plotOrder(1:(nTrigs-1)/2)),2);
wAmpsAtt(2,:) = mean(wAmps(:,plotOrder(end-(nTrigs-1)/2):end-1),2);
attNames = {'attT1','attT2'};

% PA means
for iTrig = 1:(nTrigs-1)/2 
    wAmpsPA(iTrig,:) = mean(wAmps(:,iTrig*2-1:iTrig*2),2);
end
PANames = {'T1p-T2p','T1a-T2p','T1p-T2a','T1a-T2a'};

% store results
A.attNames = attNames;
A.PANames = PANames;
A.wBaselineWindow = wBaselineWindow;
A.wAmps = wAmps;
A.wAmpsAtt = wAmpsAtt;
A.wAmpsPA = wAmpsPA;

fH = [];
fH(1) = figure;
set(gcf,'Position',tsFigPos)
set(gca,'ColorOrder',trigColors)
hold all
plot(t, wAmps(:,plotOrder))
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
plot(t, mean(wAmps(:,plotOrder(1:(nTrigs-1)/2)),2),'color',trigBlue,'LineWidth',4)
plot(t, mean(wAmps(:,plotOrder(end-(nTrigs-1)/2):end-1),2),'color',trigRed,'LineWidth',4)
legend(trigNames(plotOrder))
xlabel('time (ms)')
ylabel('wavelet amp')
title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels)])

% condition subplots
fH(2) = figure;
set(gcf,'Position',condFigPos)
for iTrig = 1:(nTrigs-1)/2
    subplot((nTrigs-1)/2,1,iTrig)
    set(gca,'ColorOrder',[trigBlue; trigRed])
    hold all
    plot(t, wAmps(:,iTrig*2-1:iTrig*2))
    legend(trigNames{iTrig*2-1:iTrig*2})
    ylim([-1 2.5])
    for iEv = 1:numel(eventTimes)
        vline(eventTimes(iEv),'k');
    end
    if iTrig==1
        title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels)])
    end
end
xlabel('time (ms)')
ylabel('wavelet amp')

% present/absent
fH(3) = figure;
set(gcf,'Position',tsFigPos)
hold on
for iTrig = 1:(nTrigs-1)/2 
    p1 = plot(t, mean(wAmps(:,iTrig*2-1:iTrig*2),2));
    set(p1, 'Color', trigColorsPA4(iTrig,:), 'LineWidth', 1.5)
end
ylim([-1 2.5])
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
legend('T1p-T2p','T1a-T2p','T1p-T2a','T1a-T2a')
xlabel('time (ms)')
ylabel('wavelet amp')
title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels)])

if saveFigs
    figPrefix = ['plot_ch' sprintf('%d_', channels) sprintf('%dHz', ssvefFreq)];
    rd_saveAllFigs(fH, {'waveletTrialAve','waveletTrialAveByCond','waveletTrialAvePA'}, figPrefix, figDir)
end

%% Hilbert on average across trials
Fbp = ssvefFreq + [-1.6 1.6];
hAmps = [];
for iTrig = 1:nTrigs
    data = trigMean(:,:,iTrig)'; % channels by samples
    dataF = ft_preproc_bandpassfilter(data,Fs,Fbp);
    dataFH = abs(hilbert(mean(dataF,1))); % average bandpassed time series across channels
    hAmps(:,iTrig) = dataFH; 
end

% attT1T2 means
hAmpsAtt(1,:) = mean(hAmps(:,plotOrder(1:(nTrigs-1)/2)),2);
hAmpsAtt(2,:) = mean(hAmps(:,plotOrder(end-(nTrigs-1)/2):end-1),2);

% PA means
for iTrig = 1:(nTrigs-1)/2 
    hAmpsPA(iTrig,:) = mean(hAmps(:,iTrig*2-1:iTrig*2),2);
end

% store results
A.hFbp = Fbp;
A.hAmps = hAmps;
A.hAmpsAtt = hAmpsAtt;
A.hAmpsPA = hAmpsPA;

fH = [];
fH(1) = figure;
set(gcf,'Position',tsFigPos)
set(gca,'ColorOrder',trigColors)
hold all
plot(t, hAmps(:,plotOrder))
legend(trigNames(plotOrder))
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
plot(t, mean(hAmps(:,plotOrder(1:(nTrigs-1)/2)),2),'color',trigBlue,'LineWidth',4)
plot(t, mean(hAmps(:,plotOrder(end-(nTrigs-1)/2):end-1),2),'color',trigRed,'LineWidth',4)
legend(trigNames(plotOrder))
xlabel('time (ms)')
ylabel('Hilbert amp')
title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels)])

% condition subplots
fH(2) = figure;
set(gcf,'Position',condFigPos)
for iTrig = 1:(nTrigs-1)/2
    subplot((nTrigs-1)/2,1,iTrig)
    set(gca,'ColorOrder',[trigBlue; trigRed])
    hold all
    plot(t, hAmps(:,iTrig*2-1:iTrig*2))
    legend(trigNames{iTrig*2-1:iTrig*2})
    for iEv = 1:numel(eventTimes)
        vline(eventTimes(iEv),'k');
    end
%     ylim([-1 2.5])
    if iTrig==1
        title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels)])
    end
end
xlabel('time (ms)')
ylabel('Hilbert amp')

% present/absent
fH(3) = figure;
set(gcf,'Position',tsFigPos)
hold on
for iTrig = 1:(nTrigs-1)/2 
    p1 = plot(t, mean(hAmps(:,iTrig*2-1:iTrig*2),2));
    set(p1, 'Color', trigColorsPA4(iTrig,:), 'LineWidth', 1.5)
end
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
legend('T1p-T2p','T1a-T2p','T1p-T2a','T1a-T2a')
xlabel('time (ms)')
ylabel('Hilbert amp')
title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels)])

% attend T1/T2 with condition error bars
mean1 = mean(hAmps(:,plotOrder(1:(nTrigs-1)/2)),2);
ste1 = std(hAmps(:,plotOrder(1:(nTrigs-1)/2)),0,2)./(sqrt((nTrigs-1)/2));
mean2 = mean(hAmps(:,plotOrder(end-(nTrigs-1)/2):end-1),2);
ste2 = std(hAmps(:,plotOrder(end-(nTrigs-1)/2):end-1),0,2)./(sqrt((nTrigs-1)/2));
fH(4) = figure;
set(gcf,'Position',tsFigPos)
hold on
shadedErrorBar(t, mean1, ste1, {'color',trigBlue,'LineWidth',4}, 1)
shadedErrorBar(t, mean2, ste2, {'color',trigRed,'LineWidth',4}, 1)
plot(t, hAmps(:,end), 'k')
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
legend('attend T1','attend T2')
xlabel('time (ms)')
ylabel('Hilbert amp')
title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels)])

if saveFigs
    figPrefix = ['plot_ch' sprintf('%d_', channels) sprintf('%dHz', ssvefFreq)];
    rd_saveAllFigs(fH, {'hilbertTrialAve','hilbertTrialAveByCond','hilbertTrialAvePA','hilbertTrialAveAttT1T2Error'}, figPrefix, figDir)
end

%% Time-frequency
taper          = 'hanning';
foi            = 1:50;
t_ftimwin      = 10 ./ foi;
toi            = tstart/1000:0.01:tstop/1000;
tfAmps = [];
for iTrig = 1:nTrigs
    data = trigMean(:,:,iTrig)'; % channels by samples
    [spectrum,ntaper,freqoi,timeoi] = ft_specest_mtmconvol(data, t/1000, ...
        'timeoi', toi, 'freqoi', foi, 'timwin', t_ftimwin, ...
        'taper', taper, 'dimord', 'chan_time_freqtap');
    specAmp = squeeze(mean(abs(spectrum),1)); % mean across channels
    tfAmps(:,:,iTrig) = specAmp';
end

tfAmpsAtt(:,:,1) = nanmean(tfAmps(:,:,plotOrder(1:(nTrigs-1)/2)),3);
tfAmpsAtt(:,:,2) = nanmean(tfAmps(:,:,plotOrder((nTrigs-1)/2+1:end-1)),3);

for iTrig = 1:(nTrigs-1)/2 
    tfAmpsPA(:,:,iTrig) = mean(tfAmps(:,:,iTrig*2-1:iTrig*2),3);
end
t1PADiff = mean(tfAmpsPA(:,:,[1 3]),3)-mean(tfAmpsPA(:,:,[2 4]),3);
t2PADiff = mean(tfAmpsPA(:,:,[1 2]),3)-mean(tfAmpsPA(:,:,[3 4]),3);

% store results
A.tfTaper = taper;
A.tfFoi = foi;
A.tfTFTimwin = t_ftimwin;
A.tfToi = toi;
A.tfAmps = tfAmps;
A.tfAmpsAtt = tfAmpsAtt;
A.tfAmpsPA = tfAmpsPA;
A.tfPADiff(:,:,1) = t1PADiff;
A.tfPADiff(:,:,2) = t2PADiff;

% figures
ytick = 10:10:numel(foi);
xtick = 51:50:numel(toi);
clims = [0 30];
diffClims = [-10 10];
hack = plotOrder;
hack(hack>4) = hack(hack>4)+1;

fH = [];
fH(1) = figure;
set(gcf,'Position',tf9FigPos)
for iTrig = 1:nTrigs
    subplot(2,5,hack(iTrig))
    imagesc(tfAmps(:,:,iTrig),clims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    if iTrig==nTrigs
        xlabel('time (ms)')
        ylabel('frequency (Hz)')
    end
    title(trigNames{iTrig})
end
rd_supertitle(['channel' sprintf(' %d', channels)]);
rd_raiseAxis(gca);

fH(2) = figure;
set(gcf,'Position',tf3FigPos)
attNames = {'attT1','attT2'};
for iAtt = 1:size(tfAmpsAtt,3)
    subplot(1,3,iAtt)
    imagesc(tfAmpsAtt(:,:,iAtt),clims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    xlabel('time (ms)')
    ylabel('frequency (Hz)')
    title(attNames{iAtt})
end
subplot(1,3,3)
imagesc(tfAmpsAtt(:,:,2)-tfAmpsAtt(:,:,1),diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (ms)')
ylabel('frequency (Hz)')
title('attT2 - attT1')
rd_supertitle(['channel' sprintf(' %d', channels)]);
rd_raiseAxis(gca);

fH(3) = figure;
set(gcf,'Position',tf9FigPos)
paNames = {'T1p-T2p','T1a-T2p','T1p-T2a','T1a-T2a'};
for iPA = 1:size(tfAmpsPA,3)
    subplot(2,4,iPA)
    imagesc(tfAmpsPA(:,:,iPA),clims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    xlabel('time (ms)')
    ylabel('frequency (Hz)')
    title(paNames{iPA})
end
subplot(2,4,5)
imagesc(t1PADiff,diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (ms)')
ylabel('frequency (Hz)')
title('T1 P-A')
subplot(2,4,6)
imagesc(t2PADiff,diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (ms)')
ylabel('frequency (Hz)')
title('T2 P-A')
subplot(2,4,7)
imagesc(t2PADiff - t1PADiff,diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (ms)')
ylabel('frequency (Hz)')
title('T2 vs. T1 P-A')
rd_supertitle(['channel' sprintf(' %d', channels)]);
rd_raiseAxis(gca);

if saveFigs
    if numel(channels)==1
        figPrefix = sprintf('im_ch%d', channels);
    else
        figPrefix = ['im_ch' sprintf('%d_', channels(1:end-1)) sprintf('%d', channels(end))];
    end
    rd_saveAllFigs(fH, {'timeFreqByCond','timeFreqAtt','timeFreqPA'}, figPrefix, figDir)
end

%% Time-frequency - single trials
taper          = 'hanning';
foi            = 1:50;
t_ftimwin      = 10 ./ foi;
toi            = tstart/1000:0.01:tstop/1000;
tfSingleAmps0 = [];
for iCh = 1:numel(channels)
    channel = channels(iCh);
    for iTrig = 1:nTrigs-1
        [iCue,iT1,iT2] = rd_indToFactorialInd(iTrig,[2,2,2]);
        data = squeeze(condData(:,channel,:,iCue,iT1,iT2))'; % trials by samples
        [spectrum,ntaper,freqoi,timeoi] = ft_specest_mtmconvol(data, t/1000, ...
            'timeoi', toi, 'freqoi', foi, 'timwin', t_ftimwin, ...
            'taper', taper, 'dimord', 'chan_time_freqtap');
        specAmp = squeeze(nanmean(abs(spectrum),1)); % mean across trials
        tfSingleAmps0(:,:,iTrig,iCh) = specAmp';
    end
end

% blank
for iCh = 1:numel(channels)
    channel = channels(iCh);
    data = squeeze(blankData(:,channel,:))';
    [spectrum,ntaper,freqoi,timeoi] = ft_specest_mtmconvol(data, t/1000, ...
        'timeoi', toi, 'freqoi', foi, 'timwin', t_ftimwin, ...
        'taper', taper, 'dimord', 'chan_time_freqtap');
    specAmp = squeeze(nanmean(abs(spectrum),1)); % mean across trials
    tfSingleAmps0(:,:,nTrigs,iCh) = specAmp';
end

% mean across channels
tfSingleAmps = mean(tfSingleAmps0, 4);

tfSingleAmpsAtt(:,:,1) = nanmean(tfSingleAmps(:,:,plotOrder(1:(nTrigs-1)/2)),3);
tfSingleAmpsAtt(:,:,2) = nanmean(tfSingleAmps(:,:,plotOrder((nTrigs-1)/2+1:end-1)),3);

for iTrig = 1:(nTrigs-1)/2 
    tfSingleAmpsPA(:,:,iTrig) = mean(tfSingleAmps(:,:,iTrig*2-1:iTrig*2),3);
end
t1SinglePADiff = mean(tfSingleAmpsPA(:,:,[1 3]),3)-mean(tfSingleAmpsPA(:,:,[2 4]),3);
t2SinglePADiff = mean(tfSingleAmpsPA(:,:,[1 2]),3)-mean(tfSingleAmpsPA(:,:,[3 4]),3);

% store results
A.stfTaper = taper;
A.stfFoi = foi;
A.stfTFTimwin = t_ftimwin;
A.stfToi = toi;
A.stfAmps = tfSingleAmps;
A.stfAmpsAtt = tfSingleAmpsAtt;
A.stfAmpsPA = tfSingleAmpsPA;
A.stfPADiff(:,:,1) = t1SinglePADiff;
A.stfPADiff(:,:,2) = t2SinglePADiff;

% figures
ytick = 10:10:numel(foi);
xtick = 51:50:numel(toi);
clims = [0 70];
diffClims = [-10 10];
hack = plotOrder;
hack(hack>4) = hack(hack>4)+1;

fH = [];
fH(1) = figure;
set(gcf,'Position',tf9FigPos)
for iTrig = 1:nTrigs
    subplot(2,5,hack(iTrig))
    imagesc(tfSingleAmps(:,:,iTrig),clims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    if iTrig==nTrigs
        xlabel('time (ms)')
        ylabel('frequency (Hz)')
    end
    title(trigNames{iTrig})
end
rd_supertitle(['channel' sprintf(' %d', channels)]);
rd_raiseAxis(gca);

fH(2) = figure;
set(gcf,'Position',tf3FigPos)
attNames = {'attT1','attT2'};
for iAtt = 1:size(tfSingleAmpsAtt,3)
    subplot(1,3,iAtt)
    imagesc(tfSingleAmpsAtt(:,:,iAtt),clims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    xlabel('time (ms)')
    ylabel('frequency (Hz)')
    title(attNames{iAtt})
end
subplot(1,3,3)
imagesc(tfSingleAmpsAtt(:,:,2)-tfSingleAmpsAtt(:,:,1),diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (ms)')
ylabel('frequency (Hz)')
title('attT2 - attT1')
rd_supertitle(['channel' sprintf(' %d', channels)]);
rd_raiseAxis(gca);

fH(3) = figure;
set(gcf,'Position',tf9FigPos)
paNames = {'T1p-T2p','T1a-T2p','T1p-T2a','T1a-T2a'};
for iPA = 1:size(tfSingleAmpsPA,3)
    subplot(2,4,iPA)
    imagesc(tfSingleAmpsPA(:,:,iPA),clims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    xlabel('time (ms)')
    ylabel('frequency (Hz)')
    title(paNames{iPA})
end
subplot(2,4,5)
imagesc(t1SinglePADiff,diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (ms)')
ylabel('frequency (Hz)')
title('T1 P-A')
subplot(2,4,6)
imagesc(t2SinglePADiff,diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (ms)')
ylabel('frequency (Hz)')
title('T2 P-A')
subplot(2,4,7)
imagesc(t2SinglePADiff - t1SinglePADiff,diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (ms)')
ylabel('frequency (Hz)')
title('T2 vs. T1 P-A')
rd_supertitle(['channel' sprintf(' %d', channels)]);
rd_raiseAxis(gca);

if saveFigs
    if numel(channels)==1
        figPrefix = sprintf('im_ch%d', channels);
    else
        figPrefix = ['im_ch' sprintf('%d_', channels(1:end-1)) sprintf('%d', channels(end))];
    end
    rd_saveAllFigs(fH, {'timeFreqSingleByCond','timeFreqSingleAtt','timeFreqSinglePA'}, figPrefix, figDir)
end

%% save analysis
if saveAnalysis
    save(analysisFileName, 'A')
end

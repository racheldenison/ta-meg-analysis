% simple_plots

% load analysis_singleTrials

%% setup
t = A.t;
eventTimes = A.eventTimes;
trigMean = A.trigMean;
nTrigs = numel(A.trigNames);
channels = A.channels;
chw = A.chw;
Fs = A.Fs;

stfFoi = A.stfFoi;
stfToi = A.stfToi;
stfAmpsAtt = A.stfAmpsAtt;
attNames = A.attNames;

%% plot setup
plotOrder = [1 5 3 7 2 6 4 8 9];
extendedMap = flipud(lbmap(nTrigs-1+4,'RedBlue'));
selectedMap = extendedMap([1:(nTrigs-1)/2 (end-(nTrigs-1)/2)+1:end],:);
trigColors = [selectedMap; 0 0 0];
trigBlue = mean(selectedMap(1:(nTrigs-1)/2,:));
trigRed = mean(selectedMap((end-(nTrigs-1)/2)+1:end,:));

tsFigPos = [0 500 1250 375];

set(0,'defaultLineLineWidth',1)

%% high-pass filtered time series
samplingInterval = 1;
tau = 100;
filtTau = samplingInterval/tau;
chIdx = 1;

yAll = [];
for iCond = 1:8
    y = squeeze(trigMean(:,chIdx,:,iCond));
    idx = isnan(y(1,:));
    y(:,idx) = [];
    yfilt = filter([1-filtTau filtTau-1],[1 filtTau-1], y);
    yMean(:,iCond) = mean(yfilt,2);
%     yMean(:,iCond) = mean(y,2);
    yAll = [yAll yfilt];
end

figure('Position',[200 50 1000 250])
plot(t,mean(yMean,2))
xlim([0 2300])
ylim([-300 300])
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
xlabel('time (ms)')
ylabel('amplitude')

figure
for i = 1:numel(A.PANames)
    subplot(4, 1, i)
    vals = yMean(:,i*2-1:i*2);
    plot(t, mean(vals,2))
    xlim([0 2300])
    for iEv = 1:numel(eventTimes)
        vline(eventTimes(iEv),'k');
    end
    title(A.PANames{i})
end

%% fft of mean time series
ym = mean(yMean,2);

% only go from cue to post-cue
% tidx1 = find(t==eventTimes(2));
% tidx2 = find(t==eventTimes(5))-1;
tidx1 = find(t==1000);
tidx2 = find(t==1500)-1;
nfft = numel(tidx1:tidx2);
Y = fft(ym(tidx1:tidx2),nfft)/nfft; % Scale by number of samples
f = Fs/2*linspace(0,1,nfft/2+1); % Fs/2 is the maximum frequency that can be measured
amps = 2*abs(Y(1:nfft/2+1)); % Multiply by 2 since only half the energy is in the positive half of the spectrum?

figure('Position',[700 150 520 220])
plot(f, amps)
xlim([0 50])
ylim([0 60])
xlabel('frequency (Hz)')
ylabel('amplitude')

%% Wavelet and ITPC
ssvefFreq = 20;

switch ssvefFreq
    case 11
        width = 4;
    case 15
        width = 6;
    case {20, 25}
        width = 8;
    case 30
        width = 12; % 12 for 30 Hz, 16 for 40 Hz gives 127 ms duration, 5 Hz bandwidth
    case 40
        width = 16;
    otherwise
        error('ssvefFreq not recognized')
end
wBaselineWindow = NaN;
% wBaselineWindow = [-500 0]; % [-300 -200];
% wBaselineWindowIdx = find(t==wBaselineWindow(1)):find(t==wBaselineWindow(2));

% only frequency of interest
wAmps0 = [];
wITPC0 = [];
wSpec0 = [];
foi = ssvefFreq;
for iTrig = 1:nTrigs
    for iCh = 1:numel(channels)
        data = squeeze(trigMean(:,iCh,:,iTrig))'; % trials by samples
        [spectrum,freqoi,timeoi] = ft_specest_wavelet(data, t/1000, 'freqoi', foi, 'width', width);
        spec = squeeze(spectrum);
        specAmp = abs(squeeze(spectrum));
        itpc = squeeze(abs(nanmean(exp(1i*angle(spectrum)),1))); % mean across trials
        
        if all(size(specAmp)>1) % if two-dimensional
            wAmp = specAmp;
            spec = spec;
        else
            wAmp = specAmp';
            spec = spec';
        end
        %     wAmpNorm = wAmp./nanmean(nanmean(wAmp(:,wBaselineWindowIdx)))-1;
        %     wAmps0(:,:,iTrig) = wAmpNorm';
        wAmps0(:,iCh,:,iTrig) = wAmp';
        wITPC0(:,iCh,iTrig) = itpc;
        wSpec0(:,:,iTrig,iCh) = spec';
    end
end
wAmps = squeeze(rd_wmean(wAmps0,chw,2)); % mean across channels
wITPC = squeeze(rd_wmean(wITPC0,chw,2));

% attT1T2 combined
att1 = []; att2 = [];
att1Spec = []; att2Spec = [];
conds1 = plotOrder(1:(nTrigs-1)/2);
conds2 = plotOrder((nTrigs-1)/2+1:nTrigs-1);
for i=1:4
    att1 = cat(2, att1, wAmps(:,:,conds1(i)));
    att2 = cat(2, att2, wAmps(:,:,conds2(i)));
    att1Spec = cat(2, att1Spec, wSpec0(:,:,conds1(i),:)); 
    att2Spec = cat(2, att2Spec, wSpec0(:,:,conds2(i),:)); 
end
wAmpsAtt(:,:,1) = att1;
wAmpsAtt(:,:,2) = att2;
wSpecAtt(:,:,1,:) = att1Spec; % time x trials x att cond x channels
wSpecAtt(:,:,2,:) = att2Spec;
attNames = {'attT1','attT2'};

% Method 2: compute itpc from all trials in condition
% attT1T2 combined
for iAtt = 1:numel(attNames)
    for iCh = 1:numel(channels)
        spectrum = wSpecAtt(:,:,iAtt,iCh)';
        itpc = squeeze(abs(nanmean(exp(1i*angle(spectrum)),1))); % mean across trials
        
        wITPCAtt0(:,iCh,iAtt) = itpc;
    end
end
wITPCAtt = squeeze(rd_wmean(wITPCAtt0,chw,2)); % mean across channels

fH = figure;
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
ylabel('wavelet amp')
title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels)])

fH = figure;
set(gcf,'Position',tsFigPos)
hold on
plot(t, wITPCAtt(:,1),'color',trigBlue,'LineWidth',2)
plot(t, wITPCAtt(:,2),'color',trigRed,'LineWidth',2)
legend(attNames)
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
xlabel('time (ms)')
ylabel('wavelet itpc')
title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels)])

%% Plot spectrum at times of interest, e.g., target-evoked peaks
tois = [-100 900 1185 1485];
nToi = numel(tois);

stfAmpsAttPeaks = [];
for iToi = 1:nToi
    toiIdx = find(abs(stfToi-tois(iToi)/1000)==min(abs(stfToi-tois(iToi)/1000)));
    
    stfAmpsAttPeaks(:,:,iToi) = squeeze(stfAmpsAtt(:,toiIdx(1),:));
end
    
ylims = [-.4 .4];
figure('Position',[350 450 600 850])
for iToi = 1:nToi
    subplot(nToi,1,iToi)
    hold on
%     plot(stfFoi([1 end]),[0 0],'--k')
    plot(stfFoi, stfAmpsAttPeaks(:,:,iToi))
    ylim(ylims)
    title(sprintf('%d ms',tois(iToi)))
    if iToi==1
        legend(attNames)
    elseif iToi==nToi
        xlabel('frequency (Hz)')
        ylabel('amplitude, normalized by trial mean')
    end
end
    
    
    
    
    
    
    
    

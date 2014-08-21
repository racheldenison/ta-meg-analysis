% rd_runERFAnalysis.m

%% Setup
filename = '/Volumes/RACHO/Data/NYU/R0890_20140806/R0890_TAPilot_8.06.14/R0890_TAPilot_8.06.14.sqd';
% filename = '/Volumes/RACHO/Data/NYU/R0817_TAPilot_8.20.14/R0817_TAPilot_8.20.14.sqd';
% trigChan = 160:167;
trigChan = 164:165; % targets
megChannels = 0:156;
channelSets = {0:39,40:79,80:119,120:156};
% badChannels = [10 11 115]; % R0890, 48-->49, 150-->152
badChannels = [10 11 115 49 152];
% badChannels = [115 152]; % R0817
tstart = -1000; % for targets
tstop = 2000;
t = tstart:tstop;

% trigNames = {'fastL-attL','fastL-attR','fastR-attL','fastR-attR',...
%     'targetL','targetR','blank'};
% trigNames = {'fastL-attL','fastL-attR','fastR-attL','fastR-attR','blank'};
trigNames = {'targetL','targetR'};

saveFigs = 0;

% load data header for plotting topologies
load data/data_hdr.mat

%% Get the data
trigMean = [];
triggers = [];
for iChSet = 1:numel(channelSets)
    allChannels = channelSets{iChSet};
    channels = setdiff(allChannels,badChannels);
    
    [trigM trigs Fs] =  rd_getData(filename, trigChan, channels, tstart, tstop);
    trigMean = cat(2,trigMean,trigM);
    triggers = cat(2,triggers,trigs); % COME BACK
end

nSamples = size(trigMean,1);
nChannels = size(trigMean,2);
nTrigs = size(trigMean,3);

%% Find noisy channels
varCutoff = 100;
trigVar = std(trigMean(:,:,1));
figure
hist(trigVar)
noisyChannels = trigVar>varCutoff;

%% Baseline
% baselinePeriod = -500:0;
baselinePeriod = t;
inBaseline = ismember(t,baselinePeriod);
baselineDC = squeeze(mean(mean(trigMean(inBaseline,:,:),1),3));
baselineTSeries = repmat(baselineDC,[size(trigMean,1),1,size(trigMean,3)]);

trigMean0 = trigMean;
trigMean = trigMean-baselineTSeries;

%% FFT on mean time series for each trigger type
% do the fft for each channel
nfft = 2^nextpow2(nSamples); % Next power of 2 from length of y
Y = fft(trigMean,nfft)/nSamples;
f = Fs/2*linspace(0,1,nfft/2+1);
amps = abs(Y(1:nfft/2+1,:,:));

%% Plot trial average and single-sided amplitude spectrum
% figure
for iTrig = 1:nTrigs
    figure
    % time
    subplot(2,1,1)
%     hold on
    plot(repmat(t',1,nChannels), trigMean(:,:,iTrig))
    xlabel('time (ms)')
    ylabel('amplitude')
    title(trigNames{iTrig})
    % frequency
    subplot(2,1,2)
%     hold on
    plot(repmat(f',1,nChannels), amps(:,:,iTrig))
    xlim([1 200])
    ylim([0 20])
    xlabel('Frequency (Hz)')
    ylabel('|Y(f)|')
end

%% Get some time series peaks
times = [228 475];
timeWindow = 5; % +/- this window value
for iT = 1:numel(times)
    time = times(iT);
    inTimeRange = t<time+timeWindow & t>time-timeWindow;
    peakMeans(iT,:,:) = mean(trigMean(inTimeRange,:,:),1);
end

%% Get the component peaks
% ssvepFreqs = [15 20 30 40];
% freqWindow = 0.2; % +/- this window value
% for iF = 1:numel(ssvepFreqs)
%     freq = ssvepFreqs(iF);
%     inFreqRange = f<freq+freqWindow & f>freq-freqWindow;
%     peakFreqs{iF} = f(inFreqRange);
%     peakVals{iF} = amps(inFreqRange,:,:);
%     peakMeans(iF,:,:) = mean(peakVals{iF},1);
% end

%% Plot peak freq image
% for iF = 1:numel(ssvepFreqs)
%     freq = ssvepFreqs(iF);
%     figure
%     imagesc(squeeze(peakMeans(iF,:,:)))
%     title(sprintf('frequency = %d',freq))
% end

%% Convert to 157 channels
timeToPlot = 228;
timeIdx = find(times==timeToPlot);
peakM = squeeze(peakMeans(timeIdx,:,:))';
inds = setdiff(0:156,badChannels)+1;
peakMeans157 = to157chan(peakM,inds,'zeros');

%% Convert to 157 channels
% freqToPlot = 40;
% freqIdx = find(ssvepFreqs==freqToPlot);
% peakM = squeeze(peakMeans(freqIdx,:,:))';
% inds = setdiff(0:156,badChannels)+1;
% peakMeans157 = to157chan(peakM,inds,'zeros');
% 
% fastLRContrast = [.5 .5 -.5 -.5 0];
% slowLRContrast = -fastLRContrast;
% fastAttInOutContrast = [.5 -.5 -.5 .5 0];
% slowAttInOutContrast = -fastAttInOutContrast;
% fastAttInOutStimLContrast = [1 -1 0 0 0];
% fastAttInOutStimRContrast = [0 0 -1 1 0];
% slowAttInOutStimLContrast = [0 0 1 -1 0];
% slowAttInOutStimRContrast = [-1 1 0 0 0];
% 
% if mod(freqToPlot,15)==0
%     lrContrast = slowLRContrast;
%     attInOutContrast = slowAttInOutContrast;
%     attInOutStimLContrast = slowAttInOutStimLContrast;
%     attInOutStimRContrast = slowAttInOutStimRContrast;
% elseif mod(freqToPlot,20)==0
%     lrContrast = fastLRContrast;
%     attInOutContrast = fastAttInOutContrast;
%     attInOutStimLContrast = fastAttInOutStimLContrast;
%     attInOutStimRContrast = fastAttInOutStimRContrast;
% end
% 
% peakStimLRDiff157 = (peakMeans157'*lrContrast')';
% peakAttInOutDiff157 = (peakMeans157'*attInOutContrast')';
% peakAttInOutDiffStimL157 = (peakMeans157'*attInOutStimLContrast')';
% peakAttInOutDiffStimR157 = (peakMeans157'*attInOutStimRContrast')';
% peakAttInOutDiffStimLRDiff157 = peakAttInOutDiffStimL157 - peakAttInOutDiffStimR157;

%% Plot on mesh
% all conditions separately
for iTrig = 1:nTrigs
    sensorData = peakMeans157(iTrig,:);
    figure
    fH = ssm_plotOnMesh(sensorData, trigNames{iTrig}, [], data_hdr, '2d');
%     set(gca,'CLim',[0 3])
%     set(gca,'CLim',[0 10])
end

% % left-right
% figure
% ssm_plotOnMesh(peakStimLRDiff157, 'L-R', [], data_hdr, '2d');
% % set(gca,'CLim',[-2 2])
% set(gca,'CLim',[-5 5])
% 
% attLims = [-3 3]; % [-.5 .5]
% % att in - att out
% figure
% ssm_plotOnMesh(peakAttInOutDiff157, 'in-out', [], data_hdr, '2d');
% set(gca,'CLim', attLims)
%     
% % left stim: att in - att out
% figure
% ssm_plotOnMesh(peakAttInOutDiffStimL157, 'L stim: in-out', [], data_hdr, '2d');
% set(gca,'CLim', attLims)
% 
% % right stim: att in - att out
% figure
% ssm_plotOnMesh(peakAttInOutDiffStimR157, 'R stim: in-out', [], data_hdr, '2d');
% set(gca,'CLim', attLims)
% 
% % att effect L stim - att effect R stim
% figure
% ssm_plotOnMesh(peakAttInOutDiffStimLRDiff157, 'att effect L - att effect R', [], data_hdr, '2d');
% set(gca,'CLim', attLims)
% 
% % save figs
% if saveFigs
%     figNames = [trigNames {'LRDiff','AttInOutDiff','LStimAttInOutDiff','RStimAttInOutDiff','AttEffectLRDiff'}];
%     rd_saveAllFigs([],figNames,sprintf('ssvef%dHz', freqToPlot))
% end

%% Find the channels with high SSVEP SNR
% peakSignal = mean(peakMeans(:,:,1:4),3); % freqs x channels
% peakNoise = mean(peakMeans(:,:,5),3); % freqs x channels
% peakSNR = (peakSignal./peakNoise)';
% peakSNRAllFlickers = mean(peakSNR,2);
% 
% figure
% imagesc(peakSNR)
% 
% figure
% hist(peakSNRAllFlickers)

%% Plot selected channels for all ssvep freqs, all conditions
% channelsToPlot = find(peakSNRAllFlickers>10);
% channelsToPlot = [14 15 26 1 50 39];
% for iF = 1:numel(ssvepFreqs)
%     freq = ssvepFreqs(iF);
%     figure
%     bar(squeeze(peakMeans(iF,channelsToPlot,:)))
%     set(gca,'XTickLabel',channelsToPlot)
%     legend(trigNames)
%     title(sprintf('%d Hz peak', freq))
% end

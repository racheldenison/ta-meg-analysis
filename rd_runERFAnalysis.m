% rd_runERFAnalysis.m

%% Setup
filename = '/Local/Users/denison/Data/TAPilot/MEG/R0890_20140806/R0890_TAPilot_8.06.14.sqd';
% filename = '/Volumes/RACHO/Data/NYU/R0890_20140806/R0890_TAPilot_8.06.14/R0890_TAPilot_8.06.14.sqd';
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
trigData = [];
for iChSet = 1:numel(channelSets)
    allChannels = channelSets{iChSet};
    channels = setdiff(allChannels,badChannels);
    
    [trigM, triggers, Fs, trigD, trigEvents] =  rd_getData(filename, trigChan, channels, tstart, tstop);
    trigMean = cat(2,trigMean,trigM);
    trigData = cat(2,trigData,trigD);
end

nSamples = size(trigMean,1);
nChannels = size(trigMean,2);
nTrigs = size(trigMean,3);

%% ANALYZE THE MEAN TIMESERIES FOR EACH TRIGGER TYPE
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
amps = 2*abs(Y(1:nfft/2+1,:,:));

%% Plot trial average and single-sided amplitude spectrum
for iTrig = 1:nTrigs
    figure
    % time
    subplot(2,1,1)
    plot(repmat(t',1,nChannels), trigMean(:,:,iTrig))
    xlabel('time (ms)')
    ylabel('amplitude')
    title(trigNames{iTrig})
    % frequency
    subplot(2,1,2)
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
    peakTMeans(iT,:,:) = mean(trigMean(inTimeRange,:,:),1);
end

%% Convert to 157 channels
timeToPlot = 475;
timeIdx = find(times==timeToPlot);
peakTM = squeeze(peakTMeans(timeIdx,:,:))';
inds = setdiff(0:156,badChannels)+1;
peakTMeans157 = to157chan(peakTM,inds,'zeros');

%% Plot on mesh
% all conditions separately
for iTrig = 1:nTrigs
    sensorData = peakTMeans157(iTrig,:);
    figure
    fH = ssm_plotOnMesh(sensorData, trigNames{iTrig}, [], data_hdr, '2d');
%     set(gca,'CLim',[0 3])
end

%% ANALYZE SINGLE TRIAL DATA
%% Baseline
% baselinePeriod = -500:0;
baselinePeriod = t;
inBaseline = ismember(t,baselinePeriod);
baselineDC = squeeze(mean(mean(trigData(inBaseline,:,:),1),3));
baselineTSeries = repmat(baselineDC,[size(trigData,1),1,size(trigData,3)]);

trigData0 = trigData;
trigData = trigData-baselineTSeries;

%% FFT on mean time series for each trigger type
% do the fft for each channel
nfft = 2^nextpow2(nSamples); % Next power of 2 from length of y
dataY = fft(trigData,nfft)/nSamples;
f = Fs/2*linspace(0,1,nfft/2+1);
dataAmps = 2*abs(dataY(1:nfft/2+1,:,:));

%% Plot the FFT for every trial from a sample channel
sampleChannel = 1;
figure
plot(f, squeeze(dataAmps(:,1,:)))
xlim([1 200])
ylim([0 40])
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')

%% Take the frequency average across trials for each channel, for each trigger type
for iTrig = 1:nTrigs
    trigger = triggers(iTrig);
    w = trigEvents(:,2)==trigger;
    dataAmpsMean(:,:,iTrig) = mean(dataAmps(:,:,w),3);
end

%% Plot single-sided amplitude spectrum for single trial data average
% figure
for iTrig = 1:nTrigs
    figure
    % frequency
    plot(repmat(f',1,nChannels), dataAmpsMean(:,:,iTrig))
    xlim([1 200])
    ylim([0 40])
    xlabel('Frequency (Hz)')
    ylabel('|Y(f)|')
    title(trigNames{iTrig})
end

%% Get the component peaks
ssvefFreqs = [15 20 30 40];
nonSSVEFFreqs = [25 25 45]; % random freqs in a similar range
freqs = [ssvefFreqs nonSSVEFFreqs];
freqWindow = 0.2; % +/- this window value
for iF = 1:numel(freqs)
    freq = freqs(iF);
    inFreqRange = f<freq+freqWindow & f>freq-freqWindow;
    peakFreqs{iF} = f(inFreqRange);
    peakVals{iF} = dataAmpsMean(inFreqRange,:,:);
    peakMeans(iF,:,:) = mean(peakVals{iF},1);
end

%% Find the channels with high SSVEF SNR
peakSignal = peakMeans(1:numel(ssvefFreqs),:,:); % freqs x channels x triggers
peakNoise = mean(peakMeans(numel(ssvefFreqs)+1:end,:,:),1); % 1 x channels x triggers
peakSNR = (peakSignal./repmat(peakNoise,numel(ssvefFreqs),1));
peakSNRAll = mean(mean(peakSNR,3),1);

snrThresh = prctile(peakSNRAll,90);
highSNRChannels = peakSNRAll>snrThresh;

% for each flicker frequency, across trigger types
figure
imagesc(mean(peakSNR,3)')

% across flicker frequencies
figure
hist(peakSNRAll)

%% Convert SNR to 157 channels and plot on mesh
peakSNRAll157 = to157chan(peakSNRAll,inds,'zeros');

figure
fH = ssm_plotOnMesh(peakSNRAll157, 'peak SNR all flickers', [], data_hdr, '2d');
set(gca,'CLim',[0 2])

%% Plot ERF of high SNR channels
for iTrig = 1:nTrigs
    figure
    % time
    subplot(2,1,1)
    plot(repmat(t',1,nnz(highSNRChannels)), trigMean(:,highSNRChannels,iTrig))
    xlabel('time (ms)')
    ylabel('amplitude')
    title(trigNames{iTrig})
    % frequency
    subplot(2,1,2)
    plot(repmat(f',1,nnz(highSNRChannels)), dataAmpsMean(:,highSNRChannels,iTrig))
    xlim([1 200])
    ylim([0 40])
    xlabel('Frequency (Hz)')
    ylabel('|Y(f)|')
end

% time, mean across channels
figure
plot(t',squeeze(mean(trigMean(:,highSNRChannels,:),2)))
xlabel('time (ms)')
ylabel('amplitude')
legend(trigNames)


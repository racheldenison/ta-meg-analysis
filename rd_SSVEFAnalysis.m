% rd_SSVEFAnalysis.m

%% Setup
exptDir = '/Local/Users/denison/Data/TAPilot/MEG';
sessionDir = 'R0890_20140806';
fileBase = 'R0890_TAPilot_8.06.14';
analStr = 'eti';

dataDir = sprintf('%s/%s', exptDir, sessionDir);

switch analStr
    case ''
        filename = sprintf('%s/%s.sqd', dataDir, fileBase);
        savename = sprintf('%s/mat/%s_ssvef_workspace.mat', dataDir, fileBase);
        figDir = sprintf('%s/figures/raw', dataDir);
    otherwise
        filename = sprintf('%s/%s_%s.sqd', dataDir, fileBase, analStr);
        savename = sprintf('%s/mat/%s_%s_ssvef_workspace.mat', dataDir, fileBase, analStr);
        figDir = sprintf('%s/figures/%s', dataDir, analStr);
end
if ~exist(figDir,'dir')
    mkdir(figDir)
end

% trigChan = 160:167;
trigChan = [160:163 166]; % stim/blank blocks
megChannels = 0:156;
channelSets = {0:39,40:79,80:119,120:156};
badChannels = [10 11 115]; % R0890
% badChannels = [115 152]; % R0817, also 152 looks dead
% badChannels = [];
tstart = 1000; % ms
tstop = 6500; % ms
t = tstart:tstop;

% trigNames = {'fastL-attL','fastL-attR','fastR-attL','fastR-attR',...
%     'targetL','targetR','blank'};
trigNames = {'fastL-attL','fastL-attR','fastR-attL','fastR-attR','blank'};

saveData = 1;
saveFigs = 1;

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

%% Save the data
if saveData
    save(savename);
end

%% Find noisy channels
varCutoff = 100;
trigVar = std(trigMean(:,:,1));
figure
hist(trigVar)
noisyChannels = trigVar>varCutoff;

%% FFT on mean time series for each trigger type
% do the fft for each channel
nfft = 2^nextpow2(nSamples); % Next power of 2 from length of y
Y = fft(trigMean,nfft)/nSamples; % Scale by number of samples
f = Fs/2*linspace(0,1,nfft/2+1); % Fs/2 is the maximum frequency that can be measured
amps = 2*abs(Y(1:nfft/2+1,:,:)); % Multiply by 2 since only half the energy is in the positive half of the spectrum?

%% Plot trial average and single-sided amplitude spectrum
% figure
fH = [];
for iTrig = 1:nTrigs
    fH(iTrig) = figure;
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

if saveFigs
    rd_saveAllFigs(fH, trigNames, 'plot_tsFFT', figDir);
end

%% Get the component peaks
ssvefFreqs = [15 20 30 40];
freqWindow = 0.2; % +/- this window value
for iF = 1:numel(ssvefFreqs)
    freq = ssvefFreqs(iF);
    inFreqRange = f<freq+freqWindow & f>freq-freqWindow;
    peakFreqs{iF} = f(inFreqRange);
    peakVals{iF} = amps(inFreqRange,:,:);
    peakMeans(iF,:,:) = mean(peakVals{iF},1);
end

%% Plot peak freq image
fH = [];
for iF = 1:numel(ssvefFreqs)
    freq = ssvefFreqs(iF);
    fH(iF) = figure;
    imagesc(squeeze(peakMeans(iF,:,:)))
    title(sprintf('frequency = %d',freq))
    freqNames{iF} = sprintf('peakAmp%dHz', freq);
end
if saveFigs
    rd_saveAllFigs(fH, freqNames, 'im', figDir)
end

%% Convert to 157 channels
for iF = 1:numel(ssvefFreqs)
    freqToPlot = ssvefFreqs(iF);
    % freqToPlot = 40;
    freqIdx = find(ssvefFreqs==freqToPlot);
    peakM = squeeze(peakMeans(freqIdx,:,:))';
    inds = setdiff(0:156,badChannels)+1;
    peakMeans157 = to157chan(peakM,inds,'zeros');
    
    fastLRContrast = [.5 .5 -.5 -.5 0];
    slowLRContrast = -fastLRContrast;
    fastAttInOutContrast = [.5 -.5 -.5 .5 0];
    slowAttInOutContrast = -fastAttInOutContrast;
    fastAttInOutStimLContrast = [1 -1 0 0 0];
    fastAttInOutStimRContrast = [0 0 -1 1 0];
    slowAttInOutStimLContrast = [0 0 1 -1 0];
    slowAttInOutStimRContrast = [-1 1 0 0 0];
    
    if mod(freqToPlot,15)==0
        lrContrast = slowLRContrast;
        attInOutContrast = slowAttInOutContrast;
        attInOutStimLContrast = slowAttInOutStimLContrast;
        attInOutStimRContrast = slowAttInOutStimRContrast;
    elseif mod(freqToPlot,20)==0
        lrContrast = fastLRContrast;
        attInOutContrast = fastAttInOutContrast;
        attInOutStimLContrast = fastAttInOutStimLContrast;
        attInOutStimRContrast = fastAttInOutStimRContrast;
    end
    
    peakStimLRDiff157 = (peakMeans157'*lrContrast')';
    peakAttInOutDiff157 = (peakMeans157'*attInOutContrast')';
    peakAttInOutDiffStimL157 = (peakMeans157'*attInOutStimLContrast')';
    peakAttInOutDiffStimR157 = (peakMeans157'*attInOutStimRContrast')';
    peakAttInOutDiffStimLRDiff157 = peakAttInOutDiffStimL157 - peakAttInOutDiffStimR157;
    
    %% Plot on mesh
    % all conditions separately
    fH = [];
    for iTrig = 1:nTrigs
        sensorData = peakMeans157(iTrig,:);
        figure
        fH(iTrig) = ssm_plotOnMesh(sensorData, trigNames{iTrig}, [], data_hdr, '2d');
        set(gca,'CLim',[0 5])
    end
    
    % left-right
    figure
    fH(end+1) = ssm_plotOnMesh(peakStimLRDiff157, 'L-R', [], data_hdr, '2d');
    set(gca,'CLim',[-3 3])
    
    attLims = [-1 1];
    % att in - att out
    figure
    fH(end+1) = ssm_plotOnMesh(peakAttInOutDiff157, 'in-out', [], data_hdr, '2d');
    set(gca,'CLim', attLims)
    
    % left stim: att in - att out
    figure
    fH(end+1) = ssm_plotOnMesh(peakAttInOutDiffStimL157, 'L stim: in-out', [], data_hdr, '2d');
    set(gca,'CLim', attLims)
    
    % right stim: att in - att out
    figure
    fH(end+1) = ssm_plotOnMesh(peakAttInOutDiffStimR157, 'R stim: in-out', [], data_hdr, '2d');
    set(gca,'CLim', attLims)
    
    % att effect L stim - att effect R stim
    figure
    fH(end+1) = ssm_plotOnMesh(peakAttInOutDiffStimLRDiff157, 'att effect L - att effect R', [], data_hdr, '2d');
    set(gca,'CLim', attLims)
    
    % save figs
    if saveFigs
        figNames = [trigNames {'LRDiff','AttInOutDiff','LStimAttInOutDiff','RStimAttInOutDiff','AttEffectLRDiff'}];
        figPrefix = sprintf('map_ssvef%dHz', freqToPlot);
        rd_saveAllFigs(fH,figNames,figPrefix,figDir)
    end
end

%% Find the channels with high SSVEF SNR
peakSignal = mean(peakMeans(:,:,1:4),3); % freqs x channels
peakNoise = mean(peakMeans(:,:,5),3); % freqs x channels
peakSNR = (peakSignal./peakNoise)';
peakSNRAllFlickers = mean(peakSNR,2);

figure
imagesc(peakSNR)

figure
hist(peakSNRAllFlickers)

%% Plot selected channels for all ssvef freqs, all conditions
% channelsToPlot = find(peakSNRAllFlickers>10);
channelsToPlot = [14 15 26 1 50 39];
for iF = 1:numel(ssvefFreqs)
    freq = ssvefFreqs(iF);
    figure
    bar(squeeze(peakMeans(iF,channelsToPlot,:)))
    set(gca,'XTickLabel',channelsToPlot)
    legend(trigNames)
    title(sprintf('%d Hz peak', freq))
end
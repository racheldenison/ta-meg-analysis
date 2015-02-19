% rd_SSVEFAnalysis.m

%% Setup
exptDir = '/Local/Users/denison/Data/TAPilot/MEG';
sessionDir = 'R0817_20140820';
fileBase = 'R0817_TAPilot_8.20.14';
analStr = 'eti'; % '', 'eti', etc.
excludeTrialsFt = 0;

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

switch sessionDir
    case 'R0890_20140806'
        % badChannels = [10 11 115]; % R0890, 48-->49, 150-->152
        badChannels = [10 11 115 49 152]; % R0890
    case 'R0817_20140820'
        badChannels = [115 152]; % R0817
        weightChannels = sort(unique([59 92 10 60 15 14 32 2 51 1 50 39 7 24 55 103 98 8]));
    otherwise
        error('sessionDir not found')
end
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
trigData = [];
trigMean = [];
triggers = [];
if excludeTrialsFt
    for iChSet = 1:numel(channelSets)
        allChannels = channelSets{iChSet};
        channels = setdiff(allChannels,badChannels);
        
        [trigM, triggers, Fs, trigD, trigEvents] =  rd_getData(filename, trigChan, channels, tstart, tstop);
        trigMean = cat(2,trigMean,trigM);
        trigData = cat(2,trigData,trigD);
    end
else
    for iChSet = 1:numel(channelSets)
        allChannels = channelSets{iChSet};
        channels = setdiff(allChannels,badChannels);
        
        [trigM trigs Fs] =  rd_getData(filename, trigChan, channels, tstart, tstop);
        trigMean = cat(2,trigMean,trigM);
        triggers = cat(2,triggers,trigs); % COME BACK
    end
end
nSamples = size(trigMean,1);
nChannels = size(trigMean,2);
nTrigs = size(trigMean,3);

%% Save the data
if saveData
    save(savename);
end

%%%%% things to set if starting from rd_prepare_vj
% nSamples = size(trigMean,1);
% nChannels = size(trigMean,2);
% nTrigs = size(trigMean,3);
% Fs = 1000;
% tstart = 1000; % ms
% tstop = 6500; % ms
% t = tstart:tstop;
% trigNames = {'fastL-attL','fastL-attR','fastR-attL','fastR-attR','blank'};
% saveFigs = 0;
% badChannels = [];
% load data/data_hdr.mat

%% Exclude trials manually rejected with ft
if excludeTrialsFt
    % load trials_rejected variable from ft manual rejection
    load([dataDir '/mat/trials_rejected_ssvef.mat'])
    
    includedTrials = logical(ones(size(trigData,3),1));
    includedTrials(trials_rejected) = 0;
    
    trigMean = [];
    for iTrig = 1:nTrigs
        trigger = triggers(iTrig);
        w = trigEvents(:,2)==trigger & includedTrials;
        trigMean(:,:,iTrig) = mean(trigData(:,:,w),3);
    end
    trigData = trigData(:,:,includedTrials);
    trigEvents = trigEvents(includedTrials,:);
    
    % update figDir
    figDir = [figDir '_ft'];
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
%     freqToPlot = 40;
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
        set(gca,'CLim',[0 15])
    end
    
    % left-right
    figure
    fH(end+1) = ssm_plotOnMesh(peakStimLRDiff157, 'L-R', [], data_hdr, '2d');
    set(gca,'CLim',[-10 10])
    
    attLims = [-3 3];
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
xlabel('peak SNR (mean across flicker frequencies)')
ylabel('number of channels')

if saveFigs
    rd_saveAllFigs(gcf, {'peakSNRAllFlickers'}, 'hist', figDir)
end

%% Find weights that maximize some conditions while minimizing others
if ~isempty(weightChannels)
    for iChannel = 1:numel(inds)
        inWeightSet(iChannel) = any(inds(iChannel)==weightChannels);
    end
else
    inWeightSet = true(ones(size(inds)));
end

condSets{1} = [1 2]; % fast left (collapse across attention conditions)
condSets{2} = [3 4]; % fast right
for iCondSet = 1:numel(condSets)
    condSet = condSets{iCondSet};
    peakMeansCond(:,:,iCondSet) = mean(peakMeans(:,inWeightSet,condSet),3);
end

% one ssvef frequency at a time
for iF = 1:numel(ssvefFreqs)
    A = squeeze(peakMeansCond(iF,:,1)); % condition 1
    B = squeeze(peakMeansCond(iF,:,2)); % condition 2
    
    % find leading generalized eigenvector of A and B
    % to maximize norm(Aw)/norm(Bw), aka maximize condition 1
    [V,D] = eig(B'*B,A'*A); % X'X so that matrices are square
    w1 = V(:,1);
    
    % to maximize norm(Bw)/norm(Aw), aka maximize condition 2
    [V,D] = eig(A'*A,B'*B);
    w2 = V(:,1);
    
    % check norm(Aw)/norm(Bw) - should be large for w1
    ((A*w1)'*(A*w1))/((B*w1)'*(B*w1)) % large
    ((B*w1)'*(B*w1))/((A*w1)'*(A*w1)) % small
    
    % check norm(Bw)/norm(Aw) - should be large for w2
    ((B*w2)'*(B*w2))/((A*w2)'*(A*w2)) % large
    ((A*w2)'*(A*w2))/((B*w2)'*(B*w2)) % small
    
    % store weights
    w(:,iF,1) = w1;
    w(:,iF,2) = w2;
end

% all frequencies simultaneously
A = peakMeansCond(:,:,1); % condition 1
B = peakMeansCond(:,:,2); % condition 2

% find leading generalized eigenvector of A and B
% to maximize norm(Aw)/norm(Bw), aka maximize condition 1
[V,D] = eig(B'*B,A'*A); % X'X so that matrices are square
w1 = V(:,1);

% to maximize norm(Bw)/norm(Aw), aka maximize condition 2
[V,D] = eig(A'*A,B'*B);
w2 = V(:,1);

% check norm(Aw)/norm(Bw) - should be large for w1
((A*w1)'*(A*w1))/((B*w1)'*(B*w1)); % large
((B*w1)'*(B*w1))/((A*w1)'*(A*w1)); % small

% check norm(Bw)/norm(Aw) - should be large for w2
((B*w2)'*(B*w2))/((A*w2)'*(A*w2)); % large
((A*w2)'*(A*w2))/((B*w2)'*(B*w2)); % small

% store weights
wAll(:,1) = w1;
wAll(:,2) = w2;

if saveData
    save(sprintf('%s/mat/weights_%s.mat',dataDir,datestr(now,'yyyymmdd')), 'w', 'wAll', 'ssvefFreqs', 'condSets', 'trigNames', 'weightChannels')
end

%% Plot weights on mesh
condNames = {'fastL','fastR'};
% inds = setdiff(0:156,badChannels)+1;
for iCond = 1:2
    w157(:,:,iCond) = to157chan(w(:,:,iCond)',weightChannels,'zeros');
end
wAll157 = to157chan(wAll',weightChannels,'zeros');

% one frequency at a time
fH = [];
for iF = 1:numel(ssvefFreqs)
    freqToPlot = ssvefFreqs(iF);
    for iCond = 1:2
        sensorData = w157(iF,:,iCond);
        figure
        fH(iCond) = ssm_plotOnMesh(sensorData, ...
            sprintf('%d Hz, %s', freqToPlot, condNames{iCond}), [], data_hdr, '2d');
        set(gca,'CLim',[-1.5 1.5])
    end
    if saveFigs
        figNames = condNames;
        figPrefix = sprintf('map_weights_ssvef%dHz', freqToPlot);
        rd_saveAllFigs(fH,figNames,figPrefix,figDir)
    end
end

% all frequencies simultaneously
fH = [];
for iCond = 1:2
    sensorData = wAll157(iCond,:);
    figure
    fH(iCond) = ssm_plotOnMesh(sensorData, condNames{iCond}, [], data_hdr, '2d');
    set(gca,'CLim',[-1.5 1.5])
end
if saveFigs
    figNames = condNames;
    figPrefix = sprintf('map_weights_allSSVEF', freqToPlot);
    rd_saveAllFigs(fH,figNames,figPrefix,figDir)
end

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
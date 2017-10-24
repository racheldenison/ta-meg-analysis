% rd_analyzeEmptyRoom.m

%% Setup
exptDir = '/Local/Users/denison/Data';
sessionDir = 'Empty_Room';
fileBase = '*9.27.17';
% fileBase = 'Shutdown';

dataDir = sprintf('%s/%s', exptDir, sessionDir);
figDir = sprintf('%s/%s', dataDir, 'figures');

runFiles = dir(sprintf('%s/%s*.sqd', dataDir, fileBase));
% filename = sprintf('%s/%s', dataDir, runFiles.name);

Fs = 1000;
channels = 0:156;
% channels = 157:159;
nChannels = numel(channels);

% badChannelIdx = [41 65 16 7];
badChannelIdx = [];

saveFigs = 1;

load data/data_hdr.mat

for iRun = 1:numel(runFiles)
filename = sprintf('%s/%s', dataDir, runFiles(iRun).name);

%% Preprocess
% preprocFileName = rd_MEGPreproc(filename, figDir, []);

%% Get data
data = sqdread(filename, 'Channels', channels);
nSamples = size(data,1);

%% Reject bad channels
data(:,badChannelIdx) = NaN;

%% segment into "trials"
trialDur = 5000;
nTrials = nSamples/trialDur;

trialData = zeros(trialDur,nChannels,nTrials);
for iChannel = 1:nChannels
    trialData(:,iChannel,:) = reshape(data(:,iChannel),trialDur,nTrials);
end

% trialData = trialData(:,:,1:10);
% nTrials = 10;

%% FFT on sigle trials
nfft = trialDur;
% nfft = 1000;
f = Fs/2*linspace(0,1,nfft/2+1); % Fs/2 is the maximum frequency that can be measured
Y = fft(trialData,nfft)/nfft; % Scale by number of samples
amps = 2*abs(Y(1:nfft/2+1,:,:)); % Multiply by 2 since only half the energy is in the positive half of the spectrum?

%% Plot spectrum
ampsMean = squeeze(nanmean(amps,3)); % mean across trials
ampsGrandMean = squeeze(nanmean(ampsMean,2)); % mean across trials and channels

toi = [15 30 60];
fH(1) = figure('Position',[600 700 450 650]);
subplot(2,1,1)
loglog(f,ampsGrandMean)
ylim([10^-1 10^3])
hold on
for i=1:numel(toi)
    vline(toi(i),'Color','k','LineStyle',':','LineWidth',1);
end
subplot(2,1,2)
loglog(f,ampsMean)
xlabel('frequency (Hz)')
ylabel('amplitude')
rd_supertitle2(und2space(rd_getTag(filename,'/')))

%% Plot ssvef
% ssvefFreqs = [12 15 20 30 40];
% ssvefFreqs = 200;
% ssvefFreqs = [14.6 27.8 29.4];
ssvefFreqs = 20;
% ssvefFreqs = 1:50;

for iFreq = 1:numel(ssvefFreqs)
    ssvefFreq = ssvefFreqs(iFreq);
    ssvefIdx = abs(f-ssvefFreq)==min(abs(f-ssvefFreq));
    freq = f(ssvefIdx);
%     ssvefIdx = f==ssvefFreq;
    
    ssvefAmps = squeeze(amps(ssvefIdx,:,:));
    ssvefAmpsMean = squeeze(nanmean(ssvefAmps,2));
    
    fH(2) = figure('Position',[700 875 900 500]);
%     clf
    subplot(1,2,1)
    imagesc(ssvefAmps)
    % set(gca,'CLim',[0 150])
    colorbar
    xlabel('trial')
    ylabel('channel')
    subplot(1,2,2)
    vals = to157chan(ssvefAmpsMean',channels+1,'zeros');
    ssm_plotOnMesh(vals, '', [], data_hdr, '2d');
    set(gca,'CLim',[0 50])
    colorbar
%     rd_supertitle2(sprintf('%s, %0.1f Hz', und2space(fileBase), freq))
    rd_supertitle2(sprintf('%s, %0.1f Hz', und2space(rd_getTag(filename,'/')), freq))
%     pause(1)

%     figure
%     vals = to157chan(ssvefAmpsMean',channels+1,'zeros');
%     ssm_plotOnMesh(vals./nanmean(vals), '', [], data_hdr, '2d');
%     set(gca,'CLim',[0 2])
%     colorbar
%     title(sprintf('%s, %d Hz', und2space(fileBase), freq))
end

%% Save figures
if saveFigs
    figNames = {'spectrum', sprintf('trials_topo_%dHz',freq)};
    rd_saveAllFigs(fH,figNames,rd_getTag(filename,'/'),figDir);
    
%     savefig(fH,sprintf('%s/%s.fig', figDir, rd_getTag(filename,'/')))
end

close all

end



% rd_runMEGPreproc.m

%% Setup
% desk
% filename = '/Local/Users/denison/Data/TAPilot/MEG/R0817_20140820/R0817_TAPilot_8.20.14.sqd';
filename = '/Local/Users/denison/Data/TAPilot/MEG/R0890_20140806/Runs/R0890_TAPilot_8.06.14_run1.sqd';

% racho
% filename = '/Volumes/RACHO/Data/NYU/R0890_20140806/R0890_TAPilot_8.06.14/R0890_TAPilot_8.06.14.sqd';
% filename = '/Volumes/RACHO/Data/NYU/R0817_TAPilot_8.20.14/R0817_TAPilot_8.20.14.sqd';

% remember these are with zero indexing
megChannels = 0:156;
refChannels = 157:159;
triggerChannels = 160:166;
photodiodeChannel = 191;

% badChannels = [];

Fl = 60; % line noise frequency
applyLineNoiseFilter = 0;

plotFigs = 0;
saveFigs = 0;

%% Get the MEG data
% data is time x channels
[data, info] = sqdread(filename);

Fs = info.SampleRate;
t = 0:1/Fs:size(data,1)/Fs-1/Fs;

%% Set aside data from special channels
% add 1 to adjust for zero-indexing
refData = data(:,refChannels+1);
trigData = data(:,triggerChannels+1);
pdData = data(:,photodiodeChannel+1);

%% Look at special channels
figure
subplot(3,1,1)
plot(t,refData)
title(['reference channels' num2str(refChannels)])
subplot(3,1,2)
plot(t,trigData)
legend(num2str(triggerChannels'))
title(['trigger channels' num2str(triggerChannels)])
subplot(3,1,3)
plot(t,pdData)
title(['photodiode channel' num2str(photodiodeChannel)])
xlabel('time (s)')

%% Denoise using reference channels
plotDenoiseFigs = 1;
% convert to time x trials x channels
data = permute(data,[1 3 2]);
data = meg_environmental_denoising(data, plotDenoiseFigs);
data = permute(data,[1 3 2]); % convert back

%% Line noise filter
if applyLineNoiseFilter
    % data should be channels x time
    data = ft_preproc_dftfilter(data', Fs, Fl);
    data = data';
end

%% Find bad channels
% high or low variance across the entire time series
outlierSDChannels = meg_find_bad_channels(permute(data(:,megChannels+1),[1 3 2]));

% dead or saturating channels for all or portions of the time series
deadChannels = checkForDeadChannels(filename);

% aggregate the bad channels
badChannels = unique([outlierSDChannels deadChannels]);
nBad = numel(badChannels);

% plot the time series for the bad channels
figure
for iBad = 1:nBad
    subplot(nBad,1,iBad)
    chan = badChannels(iBad);
    plot(t, data(:,chan))
    title(sprintf('channel %d', chan))
end
xlabel('time (s)')

%% Interpolate to replace bad channels


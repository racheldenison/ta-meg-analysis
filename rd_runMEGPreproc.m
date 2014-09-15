% rd_runMEGPreproc.m

%% Setup
% filename = '/Volumes/RACHO/Data/NYU/R0890_20140806/R0890_TAPilot_8.06.14/R0890_TAPilot_8.06.14.sqd';
% filename = '/Volumes/RACHO/Data/NYU/R0817_TAPilot_8.20.14/R0817_TAPilot_8.20.14.sqd';
filename = '/Local/Users/denison/Data/TAPilot/MEG/R0817_20140820/R0817_TAPilot_8.20.14.sqd';
% trigChan = 160:167;
trigChan = [160:163 166]; % stim/blank blocks
megChannels = 0:156;
channelSets = {0:39,40:79,80:119,120:156};
badChannels = [];
refChannels = 157:159;
photodiodeChannel = 191;
tstart = 1000; % ms
tstop = 6500; % ms
t = tstart:tstop;
Fl = 60; % line noise frequency

plotFigs = 0;
saveFigs = 0;

% load data header for plotting topologies
% load data/data_hdr.mat

%% Get the MEG data
data = [];
trigEvents = [];
for iChSet = 1:numel(channelSets)
    allChannels = channelSets{iChSet};
    channels = setdiff(allChannels,badChannels);
    
    [ev, d, Fs] = rd_get_trialsz(...
        filename, trigChan, channels, tstart, tstop);
    data = cat(2, data, d);
    trigEvents = cat(2, trigEvents, ev(:,2));
end

nSamples = size(data,1);
nChannels = size(data,2);
nTrials = size(data,3);

%% Get the reference channels and photodiode data
% Reference
[ev, refData] = rd_get_trialsz(...
    filename, trigChan, refChannels, tstart, tstop);

% Photodiode
[ev, pdData] = rd_get_trialsz(...
    filename, trigChan, photodiodeChannel, tstart, tstop);

%% Denoise using reference channels
% concatenate data and ref data and convert to time x trials x channels
ts = permute(cat(2,data,refData),[1 3 2]);
ts = meg_environmental_denoising(ts,1);
ts = ts(:,:,1:nChannels); % remove reference channels

%% Filter

%% Demean

%% Notch filter
% data should be channels x time
for iTrial = 1:nTrials
    dat0 = permute(squeeze(ts(:,iTrial,:)),[2 1]);
    dat = ft_preproc_dftfilter(dat0, Fs, Fl);
    ts(:,iTrial,:) = permute(dat, [2 1]);
end

%% Plots
if plotFigs
    % Plot photodiode
    figure
    plot(squeeze(mean(pdData,3)))
    
end


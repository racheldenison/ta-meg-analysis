% rd_runMEGPreproc.m

%% Setup
% filename = '/Volumes/RACHO/Data/NYU/R0890_20140806/R0890_TAPilot_8.06.14/R0890_TAPilot_8.06.14.sqd';
filename = '/Volumes/RACHO/Data/NYU/R0817_TAPilot_8.20.14/R0817_TAPilot_8.20.14.sqd';
% trigChan = 160:167;
trigChan = [160:163 166]; % stim/blank blocks
megChannels = 0:156;
channelSets = {0:39,40:79,80:119,120:156};
badChannels = [];
tstart = 1000; % ms
tstop = 6500; % ms
t = tstart:tstop;

% trigNames = {'fastL-attL','fastL-attR','fastR-attL','fastR-attR',...
%     'targetL','targetR','blank'};
trigNames = {'fastL-attL','fastL-attR','fastR-attL','fastR-attR','blank'};

saveFigs = 0;

% load data header for plotting topologies
load data/data_hdr.mat

%% Get the data
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


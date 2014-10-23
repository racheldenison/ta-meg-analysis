% rd_prepareFtVisRejectData.m
%
% convert visual rejection data from fieldtrip into something that we can
% load into our existing processing stream

%% Setup
exptDir = '/Local/Users/denison/Data/TAPilot/MEG';
sessionDir = 'R0817_20140820';
fileBase = 'R0817_TAPilot_8.20.14';
analStr = 'eti';
eventType = 'ssvef';

dataDir = sprintf('%s/%s', exptDir, sessionDir);

vjFileBase = sprintf('%s/mat/%s_%s_prep_%s_vj', dataDir, sessionDir, analStr, eventType);
saveFile = sprintf('%s/mat/%s_%s_%s_vjdatai.mat', dataDir, sessionDir, analStr, eventType);

switch analStr
    case ''
        figDir = sprintf('%s/figures/raw_ft', dataDir);
    otherwise
        figDir = sprintf('%s/figures/%s_ft', dataDir, analStr);
end
if ~exist(figDir,'dir')
    mkdir(figDir)
end

% load data header for plotting topologies
load data/data_hdr.mat

saveFigs = 1;

%% Load data and trigger info
load(sprintf('%s.mat', vjFileBase))
load(sprintf('%s_trl_info.mat', vjFileBase))

data = clean_SSVEP_data4;
trl = clean_SSVEP_data4_trigger_info;

clear clean_SSVEP_data4
clear clean_SSVEP_data4_trigger_info

%% Get Fs
Fs = data.fsample;

%% Make badChannels
% these will be used only for plotting
rejectedChannels = data.channels_rejected;
for iChan = 1:numel(rejectedChannels)
    badChannels(iChan) = str2double(rejectedChannels{iChan}(end-2:end));
end
badChannels = badChannels-1; % convert to 0-index

%% Visualize rejected channels
cfg = [];
cfg.layout = ft_prepare_layout(cfg, data_hdr);
cfg.interpolation = 'nearest';
cfg.style = 'straight';
cfg.electrodes = 'numbers';
cfg.colorbar = 'no';
cfg.colormap = colormap('gray');
% cfg.maplimits = 'maxmin';

rejected = ones(157, 1);
rejected(badChannels+1) = 0;

topoplot(cfg, rejected)

if saveFigs
    rd_saveAllFigs(gcf,{'rejectedChannels'},'map',figDir)
end

%% Interpolate to replace bad channels
cfg = [];
cfg.badchannel = data.channels_rejected;
cfg.method = 'spline';
data = ft_channelrepair(cfg, data);

% ft_channelrepair appends repaired channels to the end. get the order to
% place them in the standard order in trigData
[labelOrdered, chanOrder] = sort(data.label);

%% Make trigEvents, trigData, trigMean
% trigEvents is trigger time, trigger type
trigEvents = [trl(:,1)-trl(:,3) trl(:,4)];
triggers = unique(trigEvents(:,2));
nTrigs = numel(triggers);

% trigData is time x channels x events
for iEvent = 1:size(trigEvents,1)
    dat = data.trial{iEvent}';
    trigData(:,:,iEvent) = dat(:,chanOrder);
end

% trigMean is time x channels x trigger types
for iTrig = 1:nTrigs
    trigger = triggers(iTrig);
    isTrig = trigEvents(:,2)==trigger;
    trigMean(:,:,iTrig) = mean(trigData(:,:,isTrig),3);
end

%% Save data
save(saveFile, 'trigEvents','triggers','nTrigs','trigData','trigMean','badChannels','Fs')


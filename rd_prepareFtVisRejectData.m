% rd_prepareFtVisRejectData.m
%
% convert visual rejection data from fieldtrip into something that we can
% load into our existing processing stream

%% Setup
exptDir = '/Local/Users/denison/Data/TAPilot/MEG';
sessionDir = 'R0890_20140806';
fileBase = 'R0890_TAPilot_8.06.14';
analStr = 'eti';
eventType = 'ssvef';

dataDir = sprintf('%s/%s', exptDir, sessionDir);

vjFileBase = sprintf('%s/mat/%s_%s_prep_%s_vj', dataDir, sessionDir, analStr, eventType);
saveFile = sprintf('%s_/mat/%s_%s_%s_vjdata.mat', dataDir, sessionDir, analStr, eventType);

%% Load data and trigger info
load(sprintf('%s.mat', vjFileBase))
load(sprintf('%s_trl_info.mat', vjFileBase))

data = clean_SSVEP_data4;
trl = clean_SSVEP_data4_trigger_info;

clear clean_SSVEP_data4
clear clean_SSVEP_data4_trigger_info

%% Make trigEvents, trigData, trigMean
% trigEvents is trigger time, trigger type
trigEvents = [trl(:,1)-trl(:,3) trl(:,4)];
triggers = unique(trigEvents(:,2));
nTrigs = numel(triggers);

% trigData is time x channels x events
for iEvent = 1:size(trigEvents,1)
    trigData(:,:,iEvent) = data.trial{iEvent}';
end

% trigMean is time x channels x trigger types
for iTrig = 1:nTrigs
    trigger = triggers(iTrig);
    isTrig = trigEvents(:,2)==trigger;
    trigMean(:,:,iTrig) = mean(trigData(:,:,isTrig),3);
end

%% Make badChannels
% these will be used only for plotting
rejectedChannels = data.channels_rejected;
for iChan = 1:numel(rejectedChannels)
    badChannels(iChan) = str2double(rejectedChannels{iChan}(end-2:end));
end
badChannels = badChannels-1; % convert to 0-index

%% Save data



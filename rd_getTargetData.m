function [trigMean triggers Fs] =  rd_getTargetData(filename, trigChan, channels, tstart, tstop)


%% get data around each trigger
% data is [samples x channels x trigger events]
[trigEvents, data, Fs] = rd_get_trialsz(...
    filename, trigChan, channels, tstart, tstop);

%% update target triggers according to the block in which they appeared
targetTriggers = [16 32];
trigEvents(:,2) = addBlockCodeToTrigger(trigEvents(:,2), targetTriggers);
    
%% unique triggers
triggers = unique(trigEvents(:,2));
nTrigs = numel(triggers);

%% average data by trigger type
for iTrig = 1:nTrigs
    trigger = triggers(iTrig);
    events = trigEvents(:,2)==trigger;
%     trigData{iTrig} = data(:,:,events);
    trigMean(:,:,iTrig) = mean(data(:,:,events),3); % [samples x channels x trigger type]
end



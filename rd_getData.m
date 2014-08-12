function [trigMean triggers Fs] =  rd_getData(filename, trigChan, channels, tstart, tstop)


%% get data around each trigger
% data is [samples x channels x trigger events]
[trigEvents, data, Fs] = rd_get_trialsz(...
    filename, trigChan, channels, tstart, tstop);

% downsample data
% data = downsample(data,2);

triggers = unique(trigEvents(:,2));
nTrigs = numel(triggers);

%% average data by trigger type
for iTrig = 1:nTrigs
    trigger = triggers(iTrig);
    events = trigEvents(:,2)==trigger;
%     trigData{iTrig} = data(:,:,events);
    trigMean(:,:,iTrig) = mean(data(:,:,events),3); % [samples x channels x trigger type]
end



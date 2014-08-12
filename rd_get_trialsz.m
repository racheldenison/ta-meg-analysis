function [triggers, data, smpr] = rd_get_trialsz(filename, trig_chan, channels, pre, post)
%triggers comes from 'all_trigger.m' and has two columns, the first gives the position of each trigger in
%the recorded data in milliseconds; the second gives the trigger code for each trigger
%datastack is a 3D matrix with each 2D matrix in the stack consisting of a
%data matrix for each epoch; the data matrix has sensors in one dimension
%and time in the other
%the inputs include a .sqd file 'filename,' the set of trigger channels,
%e.g., [160:167] for Mac-run experiments, a list of MEG channels for
%analysis, a pre-trigger and post-trigger number for the data epochs (both
%could be positive, meaning an epoch entirely after the trigger, and both
%could be negative, meaning an epoch entirely before the trigger)
%
%
info = sqdread(filename,'info');
triggers = all_trigger(filename, trig_chan);
%gets the triggers and the trigger codes from all_triggers.m
%
%
info
%returns info about data acquisition, including sampling rate, gains, etc.
%
%converts pre and post from milliseconds to samples, e.g., if pre = 100 and
%sampling rate = 250Hz, then pre in samples is 25
pre = pre * (info.SampleRate / 1000);
post = post * (info.SampleRate / 1000);
smpr = info.SampleRate;
datastack = zeros((post-pre+1), length(channels), length(triggers(:,1)));
%creates data matrix of the correct size
%
%
for i = 1:length(triggers(:,1))
    smpl = triggers(i,1);
    epoch = sqdread(filename, 'Channels', channels, 'Samples', [smpl+pre smpl+post]);
    datastack(:,:,i) = epoch(:, :);
end;
%
%loads datastack from data
% data = baseline(datastack, pre, smpr); % RD: results in NaNs
data = datastack;
end

function [saturatedChannelEpochs, saturatedChannels, saturatedTrials] = rd_findSaturatedChannelEpochs(trigData)

cutoff = .10;

nChannels = size(trigData,2);
nTrials = size(trigData,3);

for iTrial = 1:nTrials
    saturatedChannelEpochs(:,iTrial) = findClips(squeeze(trigData(:,:,iTrial)));
end

nSaturatedTrialsEachChannel = sum(saturatedChannelEpochs,2);
nSaturatedChannelsEachTrial = sum(saturatedChannelEpochs,1);

saturatedChannels = nSaturatedTrialsEachChannel > nTrials*cutoff;
saturatedTrials = nSaturatedChannelsEachTrial > nChannels*cutoff;

figure('position',[100 100 1000 800])
colmap = colormap('parula');
hC=subplot(2,2,2);
hY=subplot(2,2,1);
hold on
hX=subplot(2,2,4);
hold on
set(hC,'Position',[.35 .35 .55 .55]);
set(hY,'Position',[.1 .35 .2 .55])
set(hX,'Position',[.35 .1 .55 .2])
bar(hY, 1:nChannels, nSaturatedTrialsEachChannel)
plot(hY, 1:nChannels, repmat(nTrials*cutoff,1,nChannels), 'k')
view(hY,-90,90) 
set(hY,'xdir','reverse')
set(hY, 'XAxisLocation', 'top')
xlim(hY, [1 nChannels])
ylim(hY, [1 nTrials])
plot(hX, 1:nTrials, nSaturatedChannelsEachTrial, 'Color', colmap(1,:))
plot(hX, 1:nTrials, repmat(nChannels*cutoff,1,nTrials), 'k')
xlim(hX, [1 nTrials])
imagesc(saturatedChannelEpochs, 'parent', hC)
xlabel(hC, 'trial')
ylabel(hC, 'channel')
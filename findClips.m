function badlyClippingChannels = findClips(data)

clippingChannelsThreshold = 4;
fractionClippingThreshold = 0.1;
numberExtraZeros = 0;

[timeOfData,numChannels] = size(data);

diffOfChannels = diff(data);
clips = diffOfChannels==0;

newClips = zeros(size(clips,1)+2,size(clips,2));
newClips( 2:end-1, :) = clips;
diffOfClips = diff(newClips);
[dummy,clipChannelEnd] = find(diffOfClips== -1);
[dummy,clipChannelBegin] = find(diffOfClips== 1);

if clipChannelEnd ~= clipChannelBegin
    error('problem');
end;

clipStart = find(diffOfClips== 1);
clipEnd = find(diffOfClips== -1);
lengthOfClips = clipEnd-clipStart;
longClips = (lengthOfClips>clippingChannelsThreshold);

clipsStartToFilter = clipStart(longClips);
clipsEndToFilter = clipEnd(longClips);

fullweights = ones(timeOfData,numChannels);
for thisClip = 1:size(clipsStartToFilter)
    fullweights(clipsStartToFilter(thisClip):clipsEndToFilter(thisClip)) = 0;
end
shiftedWeights = ones(timeOfData+numberExtraZeros*2,numChannels,numberExtraZeros*2);
for loop = 1:numberExtraZeros*2+1
    shiftedWeights(loop:end-numberExtraZeros*2-1+loop,:,loop) = fullweights;
end
paddedWeights = min(shiftedWeights,[],3);
output = paddedWeights(numberExtraZeros+1:end-numberExtraZeros,:);
totalClips = sum(output ==0,1);
fractionOfTimeClipping = totalClips/timeOfData;

badlyClippingChannels = fractionOfTimeClipping > fractionClippingThreshold;

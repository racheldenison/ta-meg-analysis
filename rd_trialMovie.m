% rd_trialMovie.m

channels = 1:40;
nTrials = size(trigData,3);

for i = 1:nTrials
    plot(squeeze(trigData(:,channels,i)))
    title(['trial ' num2str(i)])
    pause(.2)
end

for i = 1:nTrials
    imagesc(squeeze(abs(trigDataDiff(:,channels,i))))
    title(['trial ' num2str(i)])
    pause(.2)
end
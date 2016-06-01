

foi = A.tfFoi;
toi = A.tfToi;

attAlphaData = groupData.ampsAtt(8:12,:,:,:);
attAlphaData = squeeze(mean(attAlphaData));
attAlphaMean = squeeze(mean(attAlphaData,3));
attAlphaDataDiff = squeeze(attAlphaData(:,1,:)-attAlphaData(:,2,:));

figure
plot(toi,attAlphaMean)

figure
hold on
shadedErrorBar(toi*1000,attAlphaMean(:,1),std(attAlphaDataDiff,0,2)./sqrt(16),'-b',1)
shadedErrorBar(toi*1000,attAlphaMean(:,2),std(attAlphaDataDiff,0,2)./sqrt(16),'-r',1)


figure
plot(toi,attAlphaDataDiff)

figure
hold on
shadedErrorBar(toi*1000,mean(attAlphaDataDiff,2),std(attAlphaDataDiff,0,2)./sqrt(16))
for iEv = 1:numel(eventTimes)
vline(eventTimes(iEv),'color','k','LineStyle',':');
end
plot([0 3000],[0 0],'k')
h = ttest(attAlphaDataDiff');
plot(toi*1000,h)
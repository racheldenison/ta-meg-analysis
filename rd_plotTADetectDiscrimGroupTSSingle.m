function rd_plotTADetectDiscrimGroupTSSingle(A, measure, subjects, groupData, groupMean, groupSte, saveFigs, figDir, figStr)

%% args
if nargin<7
    saveFigs = 0;
end
if nargin<9
    figStr = '';
    if saveFigs==1
        error('If youre saving figs, you should specify a figStr')
    end
end

figTitle = und2space(figStr);

% groupData, A
t = A.t;
twin = A.targetWindow;
nSubjects = numel(subjects);

%% remove baseline
bwin.tsAmpsS = [-500 0];
bwin.paauTS = [twin(1) 0];

twindow.tsAmpsS = t;
twindow.paauTS = twin(1):twin(end);

fieldNames = fields(bwin);
nFields = numel(fieldNames);

for iF = 1:nFields
    fieldName = fieldNames{iF};
    tw = twindow.(fieldName);
    bw = bwin.(fieldName);
    vals = groupData.(fieldName);
    
    bidx = find(bw(1)==tw):find(bw(2)==tw);
    baseline = mean(vals(bidx,:,:,:,:));
    
    sz = size(vals);
    b = repmat(baseline,[sz(1) ones(1,length(sz)-1)]);
    
    groupDataBS.(fieldName) = vals-b; % baselined, single trials
    groupDataB.(fieldName) = squeeze(nanmean(groupDataBS.(fieldName),2)); % mean across trials
    groupMeanB.(fieldName) = mean(groupDataB.(fieldName),length(sz)-1); % mean across subjects
    groupSteB.(fieldName) = std(groupDataB.(fieldName),0,length(sz)-1)/sqrt(nSubjects); % ste across subjects
end

%% P-A 
groupDataB.paDiff = groupDataB.paauTS(:,1:2,:,:) - groupDataB.paauTS(:,3:4,:,:);
groupMeanB.paDiff = mean(groupDataB.paDiff,4);
groupSteB.paDiff = std(groupDataB.paDiff,0,4)/sqrt(nSubjects);

%% A/U ave
groupDataB.au(:,1,:,:) = mean(groupDataB.paauTS(:,[1 3],:,:),2);
groupDataB.au(:,2,:,:) = mean(groupDataB.paauTS(:,[2 4],:,:),2);
groupMeanB.au = mean(groupDataB.au,4);
groupSteB.au = std(groupDataB.au,0,4)/sqrt(nSubjects);

% auDiff
groupDataB.auDiff = -squeeze(diff(groupDataB.au,1,2));
for iT = 1:2
    [h(:,iT) p(:,iT)] = ttest(squeeze(groupDataB.auDiff(:,iT,:))');
end
groupDataB.adt = squeeze(mean(groupDataB.auDiff,2));
[hh pp] = ttest(groupDataB.adt');

%% stats on paauTS
vals = [];
vals(:,1:4,:) = groupDataB.paauTS(:,:,1,:);
vals(:,5:8,:) = groupDataB.paauTS(:,:,2,:);
condNames = {'T1_P_att','T1_P_unatt','T1_A_att','T1_A_unatt',...
    'T2_P_att','T2_P_unatt','T2_A_att','T2_A_unatt'}; 
factorNames = {'T','PA','AU'};
nLevels = [2 2 2];

for it = 1:size(vals,1)
    data = squeeze(vals(it,:,:))'; % subjects x conds
    [fvals(it,:), pvals(it,:), rowNames] = rd_rmANOVA(data, condNames, factorNames, nLevels);
end

figure
subplot(4,1,1:3)
plot(twindow,fvals)
ylabel('F value')
legend(rowNames)
subplot(4,1,4)
plot(twindow,pvals<.05)
ylim([0 2])
xlabel('time (ms)')
ylabel('p < .05')

%% plot paauT, paDiff, au
colors = get(gca,'ColorOrder');

tw = twindow.paauTS;

% paauT
valsMean = groupMeanB.paauTS;
valsSte = groupSteB.paauTS;
figure
for iT = 1:2
    subplot(2,1,iT)
    plot(tw, valsMean(:,:,iT))
    hold on
    for iPAAU = 1:4
    	shadedErrorBar(tw, valsMean(:,iPAAU,iT), valsSte(:,iPAAU,iT), {'color',colors(iPAAU,:),'LineWidth',2}, 1)
    end
    vline(0,'k');
    xlabel('time (ms)')
    ylabel('amplitude')
    title(sprintf('T%d',iT))
end
legend('P-att','P-unatt','A-att','A-unatt')

% paDiff
ylims = [-100 100];
valsMean = groupMeanB.paDiff;
valsSte = groupSteB.paDiff;
figure
for iT = 1:2
    subplot(2,1,iT)
    plot(tw, valsMean(:,:,iT))
    hold on
    for iAU = 1:2
    	shadedErrorBar(tw, valsMean(:,iAU,iT), valsSte(:,iAU,iT), {'color',colors(iAU,:),'LineWidth',2}, 1)
    end
    ylim(ylims)
    vline(0,'k');
    xlabel('time (ms)')
    ylabel('amplitude')
    title(sprintf('T%d',iT))
end
legend('att','unatt')
rd_supertitle2('P-A')

% au
ylims = [-70 70];
valsMean = groupMeanB.au;
valsSte = groupSteB.au;
figure
for iT = 1:2
    subplot(2,1,iT)
    plot(tw, valsMean(:,:,iT))
    hold on
    for iAU = 1:2
    	shadedErrorBar(tw, valsMean(:,iAU,iT), valsSte(:,iAU,iT), {'color',colors(iAU,:),'LineWidth',2}, 1)
    end
    ylim(ylims)
    vline(0,'k');
    plot(tw,h(:,iT)*10+ylims(1),'k')
    xlabel('time (ms)')
    ylabel('amplitude')
    title(sprintf('T%d',iT))
end
legend('att','unatt')

% adt
figure
hold on
plot(tw,mean(groupDataB.adt,2),'k')
shadedErrorBar(tw,mean(groupDataB.adt,2),std(groupDataB.adt,0,2)/4)
plot(twin,[0 0],'k')
plot(tw,hh*2,'k')
vline(0,'k');
xlabel('time (ms)')
ylabel('amplitude difference (att-unatt)')
title('T1 & T2')




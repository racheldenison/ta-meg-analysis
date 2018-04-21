% rd_plotTADetectDiscrimGroupDVA.m

%% setup
subjects = {'R0817_20150504', 'R0973_20150727', 'R0974_20150728', ...
    'R0861_20150813', 'R0504_20150805', 'R0983_20150813', ...
    'R0898_20150828', 'R0436_20150904', 'R1018_20151118', ...
    'R1019_20151118','R1021_20151120','R1026_20151211', ...
    'R0852_20151211','R1027_20151216','R1028_20151216',...
    'R1029_20151222'}; % N=16

subsToInclude = [1:3 5 7:16];
nSubjects = numel(subsToInclude);

channelStr = 'topChannels10'; % 'topChannels10', 'wholebrain' 
trialSelection = 'incorrect'; %'all','correct','incorrect';
respTargetSelection = 'T2Resp'; %'','T1Resp','T2Resp'

%% get data
dva = []; L2 = [];
for iSubject = 1:nSubjects
    subject = subjects{subsToInclude(iSubject)};
    filename = sprintf('/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG/%s/mat/analysis_singleTrials_*_ebi_ft_%s_%sTrials%s_dva.mat', subject, channelStr, trialSelection, respTargetSelection);
    file = dir(filename);
    load(sprintf('/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG/%s/mat/%s', subject, file.name));
    
    dva(:,iSubject) = mean(A.dva,2);
    L2(:,iSubject) = mean(A.L2,2);
end

t = A.t;
eventTimes = A.eventTimes;

%% plot figs
figure
hold on
plot(dva)
plot(mean(dva,2),'k','LineWidth',3)
for iEv = 1:numel(eventTimes)
    vline(find(t==eventTimes(iEv)),'k');
end
xlabel('time (ms)')
ylabel('dva')
title(sprintf('%s trials, N = %d', A.trialSelection, nSubjects))

figure
hold on
plot(L2)
plot(mean(L2,2),'k','LineWidth',3)
for iEv = 1:numel(eventTimes)
    vline(find(t==eventTimes(iEv)),'k');
end
xlabel('time (ms)')
ylabel('L2 norm')
title(sprintf('%s trials, N = %d', A.trialSelection, nSubjects))

%% store group data
switch trialSelection
    case 'correct'
        dvaCorrect = dva; L2Correct = L2;
    case 'incorrect'
        dvaIncorrect = dva; L2Incorrect = L2;
end

%% compare correct and incorrect
figure
subplot(2,1,1)
hold on
plot(t, [mean(dvaCorrect,2) mean(dvaIncorrect,2)])
shadedErrorBar(t, mean(dvaCorrect,2), normSte(dvaCorrect), 'b', 1)
shadedErrorBar(t, mean(dvaIncorrect,2), normSte(dvaIncorrect), 'r', 1)
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
xlabel('time (ms)')
ylabel('dva')
legend('correct','incorrect')
title(respTargetSelection)

subplot(2,1,2)
hold on
plot(t, [mean(L2Correct,2) mean(L2Incorrect,2)])
shadedErrorBar(t, mean(L2Correct,2), normSte(L2Correct), 'b', 1)
shadedErrorBar(t, mean(L2Incorrect,2), normSte(L2Incorrect), 'r', 1)
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
xlabel('time (ms)')
ylabel('L2')
legend('correct','incorrect')
title(respTargetSelection)

%% get mean of interesting time window
dvaDiff = dvaIncorrect-dvaCorrect;
[h, p, ci, stat] = ttest(dvaDiff');

switch respTargetSelection
    case 'T1Resp'
        w = t>=1550 & t<=1750;
    case 'T2Resp'
        w = t>=2000 & t<=2250;
    otherwise
        w = [];
end
dvaDiffWinMean = mean(dvaDiff(w,:),1);

figure
hold on
plot(t,dvaDiff)
plot(t,zeros(size(t)),'--k')
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
xlabel('time (ms)')
ylabel('incorrect-correct dva')

%% load behavior
b = load('/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/Behavior/Group/groupBehav_N16.mat');

switch respTargetSelection
    case 'T1Resp'
        iT = 1;
    case 'T2Resp'
        iT = 2;
    otherwise
        iT = [];
end
acc = b.groupData.overallAcc;
acc = squeeze(acc(:,iT,:));
accDiff = acc(1,:)-acc(2,:);
accDiff = accDiff(subsToInclude);

figure
twin = t(w);
scatter(accDiff, dvaDiffWinMean)
xlabel('accuracy (correct-incorrect)')
ylabel(sprintf('dva (incorrect-correct), %d-%d ms', twin(1), twin(end)))
[r p] = corr(accDiff', dvaDiffWinMean')





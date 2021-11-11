% rd_analyzeTANoiseTA2GroupBehav.m

%% setup
expt = 'TANoise'; %'TANoise','TA2'

switch expt
    case 'TANoise'
        exptDir = pathToTANoise('Behavior');
        
        %         'R1535_20190717','R1535_20190718',...
        subjects = {'R0817_20171212','R0817_20171213',...
            'R0898_20180112','R0898_20180116',...
            'R0959_20180219','R0959_20180306',...
            'R0983_20180111','R0983_20180112',...
            'R1021_20180208','R1021_20180212',...
            'R1103_20180213','R1103_20180215',...
            'R1187_20180105','R1187_20180108',...
            'R1373_20190723','R1373_20190725',...
            'R1452_20190717','R1452_20190718',...
            'R1507_20190702','R1507_20190705',...
            }; % N=10 TANoise
        
        condNames = {'valid','invalid'};
        
    case 'TA2'
        exptDir = pathToTA2('Behavior');
        
        %     'R0890_20181121', ...
        %     'R1535_20190708', 'R1535_20190711',...
        subjects = {'R0817_20181120', 'R0817_20190625',...
            'R0898_20190723', 'R0898_20190724',...
            'R0959_20181128', 'R0959_20190703',...
            'R0983_20190722', 'R0983_20190723',...
            'R1103_20181121', 'R1103_20190710',...
            'R1187_20181119', 'R1187_20190703',...
            'R1373_20181128', 'R1373_20190708',...
            'R1452_20181119', 'R1452_20190711',...
            'R1507_20190621', 'R1507_20190627',...
            'R1547_20190729', 'R1547_20190730'}; % N=10 TA2
        
        condNames = {'valid','neutral','invalid'};
end

analysisDir = sprintf('%s/Group/mat', exptDir);
analysisName = 'gN10_behavior_workspace';
saveAnalysis = 1;

% subjects = subjects([1 2 4 5 7 8 10 12 14 16]);

startRuns = ones(1,numel(subjects));

nSubjects = numel(subjects);

%% load data
for iSubject = 1:nSubjects
    sessionDir = subjects{iSubject}; 
    behavDir = sprintf('%s/%s/analysis', exptDir, sessionDir);
    behavFile = dir(sprintf('%s/*%d*.mat', behavDir, startRuns(iSubject)));
    b = load(sprintf('%s/%s', behavDir, behavFile.name));
    behav(iSubject) = behavior(b); % update behav with more info
end

%% reanalyze data
targets = unique(behav(1).responseTarget(behav(1).responseTarget~=0 & ~isnan(behav(1).responseTarget)));
cueValidities = unique(behav(1).cueValidity(~isnan(behav(1).cueValidity)));
cueValidities = sort(cueValidities,'descend'); % 1=valid, -1=invalid

groupDataAll = [];
for iS = 1:nSubjects
    for iT = 1:numel(targets)
        wT = behav(iS).responseTarget==targets(iT);
        for iV = 1:numel(cueValidities)
            wV = behav(iS).cueValidity==cueValidities(iV);
            w = wT & wV;
            groupDataAll(iS).detectHMFC{iV,iT} = behav(iS).detectHMFC(w,:);
            groupDataAll(iS).discrimCI{iV,iT} = behav(iS).discrimCI(w,:);
            groupDataAll(iS).discrimHMFC{iV,iT} = behav(iS).discrimHMFC(w,:);
            groupDataAll(iS).acc{iV,iT} = behav(iS).acc(w,:);
            groupDataAll(iS).rt{iV,iT} = behav(iS).rt(w,:);
        end
    end
end

for iS = 1:nSubjects
    detectHMFC = groupDataAll(iS).detectHMFC;
    discrimHMFC = groupDataAll(iS).discrimHMFC;
    discrimCI = groupDataAll(iS).discrimCI;
    acc = groupDataAll(iS).acc;
    rt = groupDataAll(iS).rt;
    for iT = 1:numel(targets)
        for iV = 1:numel(cueValidities)
            presentResponse = any(discrimCI{iV,iT},2);
            groupData.discrim(iV,iT,iS) = nanmean(discrimCI{iV,iT}(:,1));
            groupData.discrim1(iV,iT,iS) = nanmean(discrimCI{iV,iT}(presentResponse,1));

            groupData.detectHit(iV,iT,iS) = nanmean(detectHMFC{iV,iT}(:,1));
            groupData.detectMiss(iV,iT,iS) = nanmean(detectHMFC{iV,iT}(:,2));
            groupData.detectFA(iV,iT,iS) = nanmean(detectHMFC{iV,iT}(:,3));
            groupData.detectCR(iV,iT,iS) = nanmean(detectHMFC{iV,iT}(:,4));

            groupData.discrimHit(iV,iT,iS) = nanmean(discrimHMFC{iV,iT}(:,1));
            groupData.discrimMiss(iV,iT,iS) = nanmean(discrimHMFC{iV,iT}(:,2));
            groupData.discrimFA(iV,iT,iS) = nanmean(discrimHMFC{iV,iT}(:,3));
            groupData.discrimCR(iV,iT,iS) = nanmean(discrimHMFC{iV,iT}(:,4));

            groupData.overallAcc(iV,iT,iS) = nanmean(acc{iV,iT});
            groupData.rt(iV,iT,iS) = nanmean(rt{iV,iT});
        end
    end
end

% calculate dprime
h = groupData.detectHit;
fa = groupData.detectFA;
h(h==1) = .99;
fa(fa==0) = .01;
groupData.detectDprime = norminv(h) - norminv(fa);
groupData.detectCrit = -.5 * (norminv(h) + norminv(fa));

h = groupData.discrimHit;
fa = groupData.discrimFA;
h(h==1) = .99;
fa(fa==0) = .01;
groupData.discrimDprime = norminv(h) - norminv(fa);
groupData.discrimCrit = -.5 * (norminv(h) + norminv(fa));

%% group summary
measures = fields(groupData);
nM = numel(measures);
for iM = 1:nM
    m = measures{iM};
    sdim = numel(size(groupData.(m))); % subject dim is always the last dim
    groupMean.(m) = mean(groupData.(m),sdim);
    groupSte.(m) = std(groupData.(m),0,sdim)./sqrt(nSubjects);
end

%% average sessions
measures = fields(groupData);
nM = numel(measures);
for iM = 1:nM
    m = measures{iM};
    groupData2.(m) = (groupData.(m)(:,:,1:2:end) + groupData.(m)(:,:,2:2:end))/2;
end

%% plot group data
measures = {'discrimDprime','overallAcc','rt'}; % 'discrim1'
nM = numel(measures);
figure('color','w')
colors = get(gca,'ColorOrder');
for iM = 1:nM
    subplot(1,nM,iM)
    m = measures{iM};
    p1 = errorbar(groupMean.(m)', groupSte.(m)','.','MarkerSize',20);
    set(p1(strcmp(condNames,'invalid')),'color',colors(2,:))
    set(p1(strcmp(condNames,'neutral')),'color',[.5 .5 .5])
    xlim([.5 2.5])
    switch m
        case {'detectDprime','discrimDprime'}
            ylim([0 2.5])
        case {'detectCrit','discrimCrit'}
            ylim([-1 1])
            hold on
            plot([.5 2.5],[0 0],'--k')
        case 'overallAcc'
            ylim([0.5 1])
        case 'rt'
            ylim([0 1.6])
        otherwise
            ylim([0 1])
    end
    title(m)
    set(gca,'XTick',[1 2])
    set(gca,'XTickLabel',{'T1','T2'})
    box off
end
legend(condNames)

%% plot individual data
indivM = {'discrimDprime'};
nM = numel(indivM);
figure('Color','w','Position',[0 0 1200 800])
for iS = 1:nSubjects
    for iM = 1:numel(indivM)
        row = ceil(iS/2);
        if mod(iS,2)
            col = iM;
        else
%             col = iM+nM+1;
            col = iM+nM;
        end
%         subplot(ceil(nSubjects/2),nM*2+1,(row-1)*(nM*2+1)+col)
        subplot(ceil(nSubjects/2),nM*2,(row-1)*(nM*2)+col)
        m = indivM{iM};
        p1 = plot(groupData.(m)(:,:,iS)');
        set(p1(strcmp(condNames,'invalid')),'color',colors(2,:))
        set(p1(strcmp(condNames,'neutral')),'color',[.5 .5 .5])
        xlim([.5 2.5])
        switch m
            case {'detectDprime','discrimDprime'}
                ylim([0 3])
            case 'critDetect'
                ylim([-1 1])
                hold on
                plot([.5 2.5],[0 0],'--k')
            case 'overallAcc'
                ylim([0.5 1])
            case 'rt'
                ylim([0 1.6])
            otherwise
                ylim([0 1])
        end
        if row==1
            title(m)
        else
%             set(gca,'XTickLabel','')
%             set(gca,'YTickLabel','')
        end
%         if col==1 || col==nM+2
            ylabel(sprintf('%s',subjects{iS}))
%         end
        set(gca,'XTick',[1 2])
        set(gca,'XTickLabel',{'T1','T2'})
        box off
    end
end
legend(condNames)


%% stats
% valid vs. invalid
for iM = 1:numel(indivM)
    m = indivM{iM};
    vi.(m) = squeeze((groupData2.(m)(1,:,:)-groupData2.(m)(2,:,:)))';
    [h, pstat.(m), ci, stat.(m)] = ttest(vi.(m));
end

%% save
if saveAnalysis
    save(sprintf('%s/%s.mat', analysisDir, analysisName))
end


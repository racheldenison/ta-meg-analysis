% rd_paauAnalysis.m

%% setup
% load paauDataT1T2, paauPresAUDiffT1T2, subjects, twin
load('/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG/Group/mat/paau_workspace_20160919.mat')

tf9FigPos = [0 250 1280 580];
nSubjects = numel(subjects);
nShuffles = 10000;

paauPresAUDiffMean = squeeze(mean(paauPresAUDiffT1T2,2));

%% permutation test: any difference between A and U following target present?
% shuffle condition labels to generate null distribution of PresAUDiff
for iShuffle = 1:nShuffles
    for iT = 1:2
        attLabel = randi(2, 1, nSubjects);
        for iS = 1:nSubjects
            shuffleData(:,1,iS,iT,iShuffle) = paauDataT1T2(:,attLabel(iS),iS,iT);
            shuffleData(:,2,iS,iT,iShuffle) = paauDataT1T2(:,3-attLabel(iS),iS,iT);
        end
    end
end
shuffleAUDiff = squeeze(diff(shuffleData,1,2)); % [time subject T1T2 shuffle]
shuffleAUDiffMean = squeeze(mean(shuffleAUDiff,2)); % [time T1T2 shuffle]

for iT = 1:2
    ci(:,:,iT) = prctile(squeeze(shuffleAUDiffMean(:,iT,:)),[2.5 97.5],2);
end

%% plot mean PresAUDiff with confidence intervals
fH = figure;
for iT = 1:2
    subplot(1,2,iT)
    hold on
    plot(twin([1 end]), [0 0], 'k:')
    plot(twin(1):twin(end), ci(:,:,iT),'b','LineWidth',2)
    shadedErrorBar(twin(1):twin(end), mean(paauPresAUDiffT1T2(:,:,iT),2), std(paauPresAUDiffT1T2(:,:,iT),0,2)/sqrt(nSubjects), {'color', 'k', 'LineWidth', 3}, 1)
    vline(0,'color','k','LineStyle',':');
    xlabel('time (ms)')
    ylabel('amplitude difference (att-unatt)')
    title(sprintf('T%d, target present trials', iT))
end

%% fit line+Gaussian to individual subject curves
fittwin = [0 600];
t = fittwin(1):fittwin(end);
fittidx = [find(twin(1):twin(end)==fittwin(1)) find(twin(1):twin(end)==fittwin(2))];
opt = optimset('MaxFunEvals',20000,'MaxIter',20000);

paramNames = {'m','b','mu','sigma','amp'}; % x
linePlusGaussian = @(x,t) x(1)*t + x(2) + normpdf(t, x(3), x(4))*x(5);
cost = @(x,y,t) sum((y - linePlusGaussian(x,t)).^2);
x0 = [0, 1, 300, 50, -20];

figure
for iS = 1:nSubjects
    for iT = 1:2
        for iAU = 1:2
            % data
            y = paauDataT1T2(fittidx(1):fittidx(2),iAU,iS,iT)';
            
            % fit
            fun = @(x)cost(x,y,t);
            [x,fval,exitflag,output] = fminsearch(fun,x0,opt);
            
            % results
            yhat = linePlusGaussian(x,t);
            
            % plot
            clf
            plot([y' yhat'])
            pause(.5)
            
            % store results
            fit.y(:,iAU,iT,iS) = y;
            fit.x(:,iAU,iT,iS) = x;
            fit.cost(iAU,iT,iS) = fval;
            fit.yhat(:,iAU,iT,iS) = yhat;
            fit.exitflag(iAU,iT,iS) = exitflag;
        end
    end
end

% calculate R2
fit.sstot = squeeze(sum((fit.y - repmat(mean(fit.y),length(t),1,1,1)).^2));
fit.R2 = 1-fit.cost./fit.sstot;
fit.paramNames = paramNames;

%% plot indiv subjects with fits
ncols = ceil(sqrt(nSubjects));
nrows = ceil(nSubjects/ncols);

for iT = 1:2
    fH = figure;
    set(gcf,'Position',tf9FigPos)
    for iS = 1:nSubjects
        subplot(nrows,ncols,iS)
        hold on
        plot(twin(1):twin(end), paauDataT1T2(:,1:2,iS,iT), 'LineWidth', 2)
        plot(t, fit.yhat(:,:,iT,iS))
        %     xlim(xlims)
        %     ylim(diffYLims)
        vline(0,'color','k','LineStyle',':');
        if iS==1
            xlabel('time (ms)')
            %         ylabel('amplitude difference (T2-T1)')
        end
        title(und2space(subjects{iS}))
    end
    legend('P-att','P-unatt')
    rd_supertitle2(sprintf('T%d', iT))
end

%% plot fit results
% fit quality
figure
subplot(3,1,1)
bar(reshape(fit.cost,4,nSubjects)')
ylabel('cost')
set(gca,'XTick',1:nSubjects)
subplot(3,1,2)
bar(reshape(fit.R2,4,nSubjects)')
ylabel('R2')
set(gca,'XTick',1:nSubjects)
subplot(3,1,3)
bar(reshape(1-fit.exitflag,4,nSubjects)')
ylabel('exitflag (1=issue)')
set(gca,'XTick',1:nSubjects)

% parameters
ylims.m = [-.005 .005];
ylims.b = [-.5 2.5];
ylims.mu = [-100 1000];
ylims.sigma = [0 300];
ylims.amp = [-500 100];
nParams = numel(paramNames);
figure
for iP = 1:nParams
    subplot(nParams,1,iP)
    bar(reshape(fit.x(iP,:,:,:),4,nSubjects)')
    ylabel(paramNames{iP})
    set(gca,'XTick',1:nSubjects)
    ylim(ylims.(paramNames{iP}))
end

%% compare parameters across conditions
% exclude subjects where the fitting is bad
badFitSubjects = squeeze(any(any(fit.exitflag==0)));
posAmpSubjects = squeeze(any(any(squeeze(fit.x(strcmp(paramNames,'amp'),:,:,:))>0)));
excludeSubjects = badFitSubjects | posAmpSubjects;
goodSubjects = find(~excludeSubjects);
% goodSubjects = 1:nSubjects;

mu = squeeze(fit.x(strcmp(paramNames,'mu'),:,:,goodSubjects));
muDiff = squeeze(diff(mu));

figure
bar(reshape(mu,4,numel(goodSubjects))')
set(gca,'XTickLabel',goodSubjects)



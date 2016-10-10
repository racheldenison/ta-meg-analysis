% rd_stfStats.m

%% load data
load('/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG/Group/mat/stf_workspace_20160928.mat')

%% setup
twin = A.stftwin;
t = A.stftwinvals;
f = A.stfFoi;
eventTimes = A.eventTimes;

cmap = flipud(lbmap(64,'RedBlue'));
% cmap = colormap;
ytick = 10:10:numel(f);
paauxtick = [11 61 111];
clims = [-.04 .04];

condNames = {'T1_P_att','T1_P_unatt','T1_A_att','T1_A_unatt',...
    'T2_P_att','T2_P_unatt','T2_A_att','T2_A_unatt'};
factorNames = {'T','PA','AU'};
nLevels = [2 2 2];

nConds = numel(condNames);
nShuffles = 5000;

%% vals is freq x time x conds x subjects
vals = [];
vals(:,:,1:4,:) = groupData.PAAUT(:,:,:,1,:);
vals(:,:,5:8,:) = groupData.PAAUT(:,:,:,2,:);

%% ANOVA for each time-frequency bin
fvals = nan(numel(f),numel(t),7);
pvals = nan(numel(f),numel(t),7);
for ifreq = 1:numel(f)
    for it = 1:numel(t)
        data = squeeze(vals(ifreq,it,:,:))'; % subjects x conds
        if any(isnan(data(:)))
            % do nothing
        else
            [fvals(ifreq,it,:), pvals(ifreq,it,:)] = rd_rmANOVA(data, condNames, factorNames, nLevels);
        end
    end
end

%% main effect and interaction differences
% empirical
paDiff = squeeze(groupData.PA(:,:,1,:)-groupData.PA(:,:,2,:)); % freq x time x subjects
auDiff = squeeze(groupData.AU(:,:,1,:)-groupData.AU(:,:,2,:));

temp(:,:,1,:) = groupData.PAAU(:,:,1,:) - groupData.PAAU(:,:,3,:); % P-A, attended
temp(:,:,2,:) = groupData.PAAU(:,:,2,:) - groupData.PAAU(:,:,4,:); % P-A, unattended
paXau = squeeze(temp(:,:,1,:) - temp(:,:,2,:)); % evoked amp, att-unatt

paDiffMean = mean(paDiff,3);
auDiffMean = mean(auDiff,3);
paXauMean = mean(paXau,3);

% set nans to zeros
paDiff(isnan(paDiff)) = 0;
auDiff(isnan(auDiff)) = 0;
paXau(isnan(paXau)) = 0;

%% max cluster sum of differences
%% stage 1: empirical & threshold 
% calculate CIs for each time-freq bin to use as threshold
for iShuffle = 1:nShuffles
    % shuffle subject means independently
    for iS = 1:nSubjects
        idx = randperm(nConds);
        valsShuffle(:,:,:,iS) = vals(:,:,idx,iS);
    end

    paDiffSh = squeeze(mean(valsShuffle(:,:,[1:2 5:6],:),3) - mean(valsShuffle(:,:,[3:4 7:8],:),3));
    auDiffSh = squeeze(mean(valsShuffle(:,:,[1 3 5 7],:),3) - mean(valsShuffle(:,:,[2 4 6 8],:),3));
    
    temp = valsShuffle(:,:,[1:2 5:6],:) - valsShuffle(:,:,[3:4 7:8],:); % Patt-Aatt:T1, Punatt-Aunatt:T1
    paXauSh = squeeze(mean(temp(:,:,[1 3],:) - temp(:,:,[2 4],:),3)); % mean across targets
    
    paDiffNull(:,:,iShuffle) = mean(paDiffSh,3); % mean across subjects
    auDiffNull(:,:,iShuffle) = mean(auDiffSh,3);
    paXauNull(:,:,iShuffle) = mean(paXauSh,3);
end

paDiffCI = prctile(paDiffNull,[2.5 97.5],3);
auDiffCI = prctile(auDiffNull,[2.5 97.5],3);
paXauCI = prctile(paXauNull,[2.5 97.5],3);

paThresh = paDiffMean>paDiffCI(:,:,2) | paDiffMean<paDiffCI(:,:,1);
auThresh = auDiffMean>auDiffCI(:,:,2) | auDiffMean<auDiffCI(:,:,1);
paXauThresh = paXauMean>paXauCI(:,:,2) | paXauMean<paXauCI(:,:,1);

% empirical cluster sum
[~, paDiffCluster] = rd_clusterSum2D(paDiffMean, paThresh);
[~, auDiffCluster] = rd_clusterSum2D(auDiffMean, auThresh);
[~, paXauCluster] = rd_clusterSum2D(paXauMean, paXauThresh);

figure
subplot(2,2,1)
imagesc(paDiffMean)
subplot(2,2,2)
imagesc(paThresh)
subplot(2,2,3)
imagesc(paDiffCI(:,:,1))
subplot(2,2,4)
imagesc(paDiffCI(:,:,2))
aH = findall(gcf,'type','axes');
for iAx = 1:numel(aH)
    axes(aH(iAx));
    rd_timeFreqPlotLabels(t,f,paauxtick,ytick,0);
    set(gca,'clim',clims)
    if iAx==2
        xlabel('time (s)')
        ylabel('frequency (Hz)')
    end
end
colormap(cmap)
rd_supertitle2('P-A')

figure
subplot(2,2,1)
imagesc(auDiffMean)
subplot(2,2,2)
imagesc(auThresh)
subplot(2,2,3)
imagesc(auDiffCI(:,:,1))
subplot(2,2,4)
imagesc(auDiffCI(:,:,2))
aH = findall(gcf,'type','axes');
for iAx = 1:numel(aH)
    axes(aH(iAx));
    rd_timeFreqPlotLabels(t,f,paauxtick,ytick,0);
    set(gca,'clim',clims)
    if iAx==2
        xlabel('time (s)')
        ylabel('frequency (Hz)')
    end
end
colormap(cmap)
rd_supertitle2('A-U')

figure
subplot(2,2,1)
imagesc(paXauMean)
subplot(2,2,2)
imagesc(paXauThresh)
subplot(2,2,3)
imagesc(paXauCI(:,:,1))
subplot(2,2,4)
imagesc(paXauCI(:,:,2))
aH = findall(gcf,'type','axes');
for iAx = 1:numel(aH)
    axes(aH(iAx));
    rd_timeFreqPlotLabels(t,f,paauxtick,ytick,0);
    set(gca,'clim',clims)
    if iAx==2
        xlabel('time (s)')
        ylabel('frequency (Hz)')
    end
end
colormap(cmap)
rd_supertitle2('paXau')

figure
hold on
plot(auDiffNull)
plot(auDiffMean,'k','LineWidth',2)
plot(auDiffCI,'LineWidth',2)
plot(5*auThresh,'k','LineWidth',2)

figure
hold on
plot(paDiffNull)
plot(paDiffMean,'k','LineWidth',2)
plot(paDiffCI,'LineWidth',2)
plot(5*paThresh,'k','LineWidth',2)

%% stage 2: max cluster sum permutation
% use CIs calculated in previous stage
for iShuffle = 1:nShuffles
    % shuffle subject means independently
    for iS = 1:nSubjects
        idx = randperm(nConds);
        valsShuffle(:,:,iS) = vals(twoi,idx,iS);
    end

    paDiffShMean = mean(squeeze(mean(valsShuffle(:,[1:2 5:6],:),2) - mean(valsShuffle(:,[3:4 7:8],:),2)),2);
    auDiffShMean = mean(squeeze(mean(valsShuffle(:,[1 3 5 7],:),2) - mean(valsShuffle(:,[2 4 6 8],:),2)),2);
    
    paThreshSh = paDiffShMean>paDiffCI(:,2) | paDiffShMean<paDiffCI(:,1);
    auThreshSh = auDiffShMean>auDiffCI(:,2) | auDiffShMean<auDiffCI(:,1);
    
    [~, paDiffClusterNull(iShuffle)] = rd_clusterSum(paDiffShMean, paThreshSh);
    [~, auDiffClusterNull(iShuffle)] = rd_clusterSum(auDiffShMean, auThreshSh);
end

paDiffClusterCI = prctile(paDiffClusterNull,[2.5 97.5],2);
auDiffClusterCI = prctile(auDiffClusterNull,[2.5 97.5],2);

%% max cluster sum of t-stat
%% stage 1: empirical & threshold
for ifreq = 1:numel(f)
    [h p ci stats] = ttest(squeeze(paDiff(ifreq,:,:))');
    paDiffTStat(ifreq,:) = stats.tstat;
    
    [h p ci stats] = ttest(squeeze(auDiff(ifreq,:,:))');
    auDiffTStat(ifreq,:) = stats.tstat;
    
    [h p ci stats] = ttest(squeeze(paXau(ifreq,:,:))');
    paXauTStat(ifreq,:) = stats.tstat;
end

% t threshold
tthresh = abs(tinv(.05/2,nSubjects-1));

% empirical cluster sum
[~, paDiffCluster] = rd_clusterSum2D(paDiffTStat, abs(paDiffTStat)>tthresh);
[~, auDiffCluster] = rd_clusterSum2D(auDiffTStat, abs(auDiffTStat)>tthresh);
[~, paXauCluster] = rd_clusterSum2D(paXauTStat, abs(paXauTStat)>tthresh);

% plot
clims = [-4 4];
figpos = [200 500 750 300];
figure('Position',figpos);
subplot(1,2,1)
imagesc(paDiffTStat)
set(gca,'clim',clims)
subplot(1,2,2)
imagesc(abs(paDiffTStat)>tthresh)
aH = findall(gcf,'type','axes');
for iAx = 1:numel(aH)
    axes(aH(iAx));
    rd_timeFreqPlotLabels(t,f,paauxtick,ytick,0);
    if iAx==1
        xlabel('time (s)')
        ylabel('frequency (Hz)')
    end
end
colormap(cmap)
rd_supertitle2('P-A')

figure('Position',figpos);
subplot(1,2,1)
imagesc(auDiffTStat)
set(gca,'clim',clims)
subplot(1,2,2)
imagesc(abs(auDiffTStat)>tthresh)
aH = findall(gcf,'type','axes');
for iAx = 1:numel(aH)
    axes(aH(iAx));
    rd_timeFreqPlotLabels(t,f,paauxtick,ytick,0);
    if iAx==1
        xlabel('time (s)')
        ylabel('frequency (Hz)')
    end
end
colormap(cmap)
rd_supertitle2('A-U')

figure('Position',figpos);
subplot(1,2,1)
imagesc(paXauTStat)
set(gca,'clim',clims)
subplot(1,2,2)
imagesc(abs(paXauTStat)>tthresh)
aH = findall(gcf,'type','axes');
for iAx = 1:numel(aH)
    axes(aH(iAx));
    rd_timeFreqPlotLabels(t,f,paauxtick,ytick,0);
    if iAx==1
        xlabel('time (s)')
        ylabel('frequency (Hz)')
    end
end
colormap(cmap)
rd_supertitle2('paXau')

%% stage 2: max cluster sum permutation
for iShuffle = 1:nShuffles
    if mod(iShuffle,100)==0
        fprintf('.')
        if mod(iShuffle,1000)==0
            fprintf('%d %s\n',iShuffle, datestr(now))
        end
    end
    % shuffle subject means independently
    for iS = 1:nSubjects
        idx = randperm(nConds);
        valsShuffle(:,:,:,iS) = vals(:,:,idx,iS);
    end

    paDiffSh = squeeze(mean(valsShuffle(:,:,[1:2 5:6],:),3) - mean(valsShuffle(:,:,[3:4 7:8],:),3));
    auDiffSh = squeeze(mean(valsShuffle(:,:,[1 3 5 7],:),3) - mean(valsShuffle(:,:,[2 4 6 8],:),3));
    
    temp = valsShuffle(:,:,[1:2 5:6],:) - valsShuffle(:,:,[3:4 7:8],:); % Patt-Aatt:T1, Punatt-Aunatt:T1
    paXauSh = squeeze(mean(temp(:,:,[1 3],:) - temp(:,:,[2 4],:),3)); % mean across targets
    
    % set nans to zeros
    paDiffSh(isnan(paDiffSh)) = 0;
    auDiffSh(isnan(auDiffSh)) = 0;
    paXauSh(isnan(paXauSh)) = 0;
    
    for ifreq = 1:numel(f)
        [h p ci stats] = ttest(squeeze(paDiffSh(ifreq,:,:))');
        paDiffShTStat(ifreq,:) = stats.tstat;
        
        [h p ci stats] = ttest(squeeze(auDiffSh(ifreq,:,:))');
        auDiffShTStat(ifreq,:) = stats.tstat;
        
        [h p ci stats] = ttest(squeeze(paXauSh(ifreq,:,:))');
        paXauShTStat(ifreq,:) = stats.tstat;
    end
    
    [~, paDiffClusterNull(iShuffle)] = rd_clusterStat(paDiffShTStat, abs(paDiffShTStat)>tthresh);
    [~, auDiffClusterNull(iShuffle)] = rd_clusterStat(auDiffShTStat, abs(auDiffShTStat)>tthresh);
    [~, paXauClusterNull(iShuffle)] = rd_clusterStat(paXauShTStat, abs(paXauShTStat)>tthresh);
end

paDiffClusterCI = prctile(paDiffClusterNull,95);
auDiffClusterCI = prctile(auDiffClusterNull,95);
paXauClusterCI = prctile(paXauClusterNull,95);

%% plot empirical clusters exceeding cluster threshold
tmap = paDiffTStat;
ci = paDiffClusterCI;
threshmap = double(abs(tmap)>tthresh);
[cstats, cmax, C] = rd_clusterStat(tmap, threshmap);

for iCluster = 1:numel(cstats)
    if abs(cstats(iCluster))>ci
        threshmap(C.PixelIdxList{iCluster}) = 2;
    end
end
figure
imagesc(threshmap)
rd_timeFreqPlotLabels(t,f,paauxtick,ytick,0);
colormap(cmap)
xlabel('time (s)')
ylabel('frequency (Hz)')


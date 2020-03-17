% rd_plotVHTimeFreq.m

%%
exptDir = pathToTANoise('MEG');
analysisDir = sprintf('%s/Group/mat', exptDir);

%% load data
measure = 'stfITPC-single';
selectionStr = 'topChannels5_allTrials';
normalizeOption = 'none';
% also set exptType and ssvefFreq in function

[groupData, groupMean, groupSte, A] = rd_plotTADetectDiscrimGroup(measure, selectionStr, normalizeOption);

%% setup
targetNames = {'T1','T2'};
attNames = {'att','unatt'};

nT = numel(targetNames);
nAtt = numel(attNames);

twin = A.wtwin;

nTrigs = size(A.ampsMean,3);
nSubjects = size(groupData.amps,4)/2; % /2 assuming 2 sessions per subject

% plot
% plotOrder = [1 5 3 7 2 6 4 8 9];

tf9SquareFigPos = [50 50 850 850];
tf6SquareFigPos = [50 50 850 530];

cmap = flipud(lbmap(64,'RedBlue'));
% cmap = colormap;

foi = A.stfFoi;
toi = A.stftwinvals;

ytick = 10:10:numel(foi);
paauxtick = [11 61 111];

switch measure
    case 'stf-single'
        clims = [-0.1 0.1];
        diffClims = [-0.07 0.07];
    case 'stfITPC-single'
        clims = [0 0.5];
        diffClims = [-0.07 0.07];
    otherwise
        error('measure not recognized')
end

switch A.exptType
    case 'TADetectDiscrim'
        PresAbsNames = {'present','absent'};
        PresAbsNamesShort = {'P','A'};
    case 'TANoise'
        PresAbsNames = {'vertical','horizontal'};
        PresAbsNamesShort = {'V','H'};
    otherwise
        error('exptType not recognized')
end

%% calculate V-H average
% combine subject sessions
vh = [];
for iS = 1:nSubjects
    vh(:,:,:,iS) = (groupData.PA(:,:,:,iS*2-1) + groupData.PA(:,:,:,iS*2))/2;
end

vhDiff = squeeze(vh(:,:,1,:) - vh(:,:,2,:));
vhDiff = vhDiff(:,find(toi==0):end,:);

%% t-test at all time-freq bins
zz = zeros(size(vhDiff));
[h, p, ci, stat] = ttest(vhDiff,zz,'dim',3);
tvals = stat.tstat;
tvals(isnan(tvals)==0);

% cluster sums
[clusterStats, maxClusterStat, C] = rd_clusterStat(abs(tvals), h);

%% permute V and H conditions for each subject to get null distribution
nP = 1000;

maxClusterStatSh = [];
for iP = 1:nP
    % flip sign of half of the subject difference maps
    x = randperm(nSubjects);
    flip = mod(x,2);
    flip(flip==0) = -1;
    
    % flip subject difference maps to shuffle
    vhDiffSh = [];
    for iS = 1:nSubjects
        vhDiffSh(:,:,iS) = vhDiff(:,:,iS).*flip(iS);
    end
    
    % t-tests
    [h, p, ci, stat] = ttest(vhDiffSh,zz,'dim',3);
    tvalsSh = stat.tstat;
    tvalsSh(isnan(tvalsSh)==0);
    
    % cluster sums
    [clusterStatsSh, maxClusterStatSh(iP)] = rd_clusterStat(abs(tvalsSh), h);
end

%% plot
figure
imagesc(mean(vhDiff,3))
rd_timeFreqPlotLabels(toi,foi,paauxtick,ytick,0);
colorbar
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('ITPC V-H')

figure
imagesc(abs(tvals))
rd_timeFreqPlotLabels(toi,foi,paauxtick,ytick,0);
colorbar
set(gca,'clim',[tthresh max(abs(tvals(:)))])
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('ITPC V-H t stat')

figure
histogram(maxClusterStatSh)
vline(maxClusterStat,'r')
xlabel('Max cluster stat')
ylabel('Frequency in shuffled data')
title('ITPC V-H t stat clusters')

% format subplots
% aH = findall(gcf,'type','axes');
% for iAx = 1:numel(aH)
%     axes(aH(iAx));
%     rd_timeFreqPlotLabels(toi,foi,paauxtick,ytick,0);
%     set(gca,'clim',diffClims)
% end
% colormap(cmap)
% rd_supertitle2('V-H')


%% save analysis
vhFileName = sprintf('%s/gN%d_VHTimeFreq_workspace.mat', analysisDir, nSubjects);

if saveAnalysis
    save(vhFileName)
end


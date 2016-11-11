function rd_plotTADetectDiscrimGroupPreWholebrain(A, measure, subjects, ...
    groupData, groupMean, groupSte, groupTStat, ...
    saveFigs, figDir, figStr)

%% setup
% load data header for plotting topologies
load data/data_hdr.mat
cfg = [];
layout = ft_prepare_layout(cfg, data_hdr);

load data/grad.mat
cfg = [];
cfg.method = 'distance';
% cfg.neighbourdist = 4.5; % default is 4
cfg.grad = grad;
cfg.feedback = 'yes';
neighbours = ft_prepare_neighbours(cfg);

load parula
cmap = flipud(lbmap(64,'RedBlue'));
clims = [-5 5];

%% plot
fieldNames = fields(groupTStat);
nFields = numel(fieldNames);

fH = [];
fH(1) = figure('Position',[360 280 450 630]);
for iF = 1:nFields
    fieldName = fieldNames{iF};
    for iT = 1:2
        subplot(nFields,2,iT+(iF-1)*2)
        str = sprintf('T%d %s',iT,fieldName);
        vals = groupTStat.(fieldName)(:,iT)';
        ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
        set(gca,'CLim',clims)
        colormap(cmap)
    end
end
rd_supertitle2('A-U t-stat')

%% cluster stats
design(1,:) = [1:nSubjects 1:nSubjects];
design(2,:) = [ones(1,nSubjects) 2*ones(1,nSubjects)];

cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'depsamplesT';
cfg.correctm = 'cluster';
cfg.clusteralpha = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan = 0; % minimum number of neighborhood channels to be included in clustering
cfg.neighbours = neighbours;
cfg.tail = 0;
cfg.clustertail = 0;
cfg.alpha = 0.025;
cfg.numrandomization = 500;
cfg.design = design;
cfg.uvar = 1; % unit variable
cfg.ivar = 2; % independent variable

figure
fieldName = 'alpha';
for iT = 1:2
    data = [];
    for iObs = 1:size(design,2)
        subject = design(1,iObs);
        cond = design(2,iObs);
        data{iObs}.label = layout.label(1:numel(A.channels));
        data{iObs}.fsample = A.Fs;
        data{iObs}.avg = squeeze(groupData.(fieldName)(:,cond,iT,subject));
        data{iObs}.time = 0;
        data{iObs}.dimord = 'chan_time';
    end
    
    stat = ft_timelockstatistics(cfg, data{:});
    
    if isfield(stat,'posclusters') && ~isempty(stat.posclusters)
        pos_cluster_pvals = [stat.posclusters(:).prob];
        pos_signif_clust = find(pos_cluster_pvals < stat.cfg.alpha);
        pos = ismember(stat.posclusterslabelmat, pos_signif_clust);
    else
        pos = [];
    end
    if isfield(stat,'negclusters') && ~isempty(stat.negclusters)
        neg_cluster_pvals = [stat.negclusters(:).prob];
        neg_signif_clust = find(neg_cluster_pvals < stat.cfg.alpha);
        neg = ismember(stat.negclusterslabelmat, neg_signif_clust);
    else
        neg = [];
    end
    
    hc = [find(pos); find(neg)];
    
    subplot(1,2,iT)
    rd_topoplot(stat.stat, layout, cmap, hc)
    set(gca,'clim',clims)
    rd_supertitle2(fieldName)
end


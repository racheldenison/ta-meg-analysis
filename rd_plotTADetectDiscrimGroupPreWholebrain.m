function rd_plotTADetectDiscrimGroupPreWholebrain(A, measure, subjects, ...
    groupData, groupMean, groupSte, groupTStat, ...
    saveFigs, figDir, figStr)

%% setup
% load data header for plotting topologies
load data/data_hdr.mat
cfg = [];
layout = ft_prepare_layout(cfg, data_hdr);

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
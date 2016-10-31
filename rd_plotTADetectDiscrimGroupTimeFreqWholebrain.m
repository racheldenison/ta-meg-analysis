function rd_plotTADetectDiscrimGroupTimeFreqWholebrain(A, measure, subjects, ...
    groupData, groupMean, groupSte, groupTStat, ...
    saveFigs, figDir, figStr)

%% set up
foi = A.stfFoi;
toi = A.stfToi;
twinvals = A.stftwinvals;
trigNames = A.trigNames;
nTrigs = numel(trigNames);

ytick = 10:10:numel(foi);
xtick = 51:50:numel(toi);
load parula
cmap = flipud(lbmap(64,'RedBlue'));

switch measure
    case 'stf-single-wb'
        clims = [-0.15 0.15];
        diffClims = [-0.1 0.1];
    case 'itpc-single-wb'
        clims = [.08 0.22];
        diffClims = [-0.02 0.02];
end

plotOrder = [1 5 3 7 2 6 4 8 9];
hack = plotOrder;
hack(hack>4) = hack(hack>4)+1;

tf9FigPos = [0 250 1280 580];
tf3FigPos = [200 475 980 330];
tf9SquareFigPos = [50 50 850 850];

load data/data_hdr.mat
cfg = [];
layout = ft_prepare_layout(cfg, data_hdr);

%% whole trial by condition
fH = [];
fH(1) = figure;
set(gcf,'Position',tf9FigPos)
for iTrig = 1:nTrigs
    subplot(2,5,hack(iTrig))
    imagesc(squeeze(nanmean(groupMean.amps(:,:,:,iTrig),1)),clims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    colormap(parula)
    if iTrig==nTrigs
        xlabel('time (s)')
        ylabel('frequency (Hz)')
    end
    title(trigNames{iTrig})
end
rd_supertitle2('amplitude, all channels mean');

if strcmp(measure,'stf-single-wb')
    fH(2) = figure;
    set(gcf,'Position',tf3FigPos)
    attNames = {'attT1','attT2'};
    for iAtt = 1:size(groupMean.ampsAtt,4)
        subplot(1,3,iAtt)
        imagesc(squeeze(nanmean(groupMean.ampsAtt(:,:,:,iAtt),1)),clims)
        title(attNames{iAtt})
        colormap(parula)
        freezeColors
    end
    subplot(1,3,3)
    imagesc(squeeze(nanmean((groupMean.ampsAtt(:,:,:,2)-groupMean.ampsAtt(:,:,:,1)),1)),diffClims)
    title('attT2 - attT1')
    aH = findall(gcf,'type','axes');
    for iAx = 1:numel(aH)
        axes(aH(iAx));
        rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
        xlabel('time (s)')
        ylabel('frequency (Hz)')
    end
    colormap(cmap)
    rd_supertitle2('amplitude, all channels mean');
    
    fH(3) = figure;
    set(gcf,'Position',tf9FigPos)
    paNames = {'T1p-T2p','T1a-T2p','T1p-T2a','T1a-T2a'};
    for iPA = 1:size(groupMean.ampsPA,4)
        subplot(2,4,iPA)
        imagesc(squeeze(nanmean(groupMean.ampsPA(:,:,:,iPA),1)),clims)
        rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
        xlabel('time (s)')
        ylabel('frequency (Hz)')
        title(paNames{iPA})
        colormap(parula)
        freezeColors
    end
    subplot(2,4,5)
    imagesc(squeeze(nanmean(groupMean.paDiff(:,:,:,1),1)),diffClims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    xlabel('time (s)')
    ylabel('frequency (Hz)')
    title('T1 P-A')
    subplot(2,4,6)
    imagesc(squeeze(nanmean(groupMean.paDiff(:,:,:,2),1)),diffClims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    xlabel('time (s)')
    ylabel('frequency (Hz)')
    title('T2 P-A')
    subplot(2,4,7)
    imagesc(squeeze(nanmean((groupMean.paDiff(:,:,:,2) - groupMean.paDiff(:,:,:,1)),1)),diffClims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    xlabel('time (s)')
    ylabel('frequency (Hz)')
    title('T2 vs. T1 P-A')
    rd_supertitle2('amplitude, all channels mean');
    colormap(cmap)
end

%% 9 squares, attended-unattended
fH(4) = figure;
set(gcf,'Position',tf9SquareFigPos)
% T1/T2 x pres/abs
subplot(3,3,1)
imagesc(squeeze(nanmean((groupMean.PAAUT(:,:,:,1,1)-groupMean.PAAUT(:,:,:,2,1)),1))) % T1-pres-att vs. unatt
ylabel('present')
title('T1')
subplot(3,3,2)
imagesc(squeeze(nanmean((groupMean.PAAUT(:,:,:,1,2)-groupMean.PAAUT(:,:,:,2,2)),1))) % T2-pres-att vs. unatt
title('T2')
subplot(3,3,4)
imagesc(squeeze(nanmean((groupMean.PAAUT(:,:,:,3,1)-groupMean.PAAUT(:,:,:,4,1)),1))) % T1-abs-att vs. unatt
ylabel('absent')
subplot(3,3,5)
imagesc(squeeze(nanmean((groupMean.PAAUT(:,:,:,3,2)-groupMean.PAAUT(:,:,:,4,2)),1))) % T2-abs-att vs. unatt
% ave(T1,T2)
subplot(3,3,3)
imagesc(squeeze(nanmean((groupMean.PAAU(:,:,:,1)-groupMean.PAAU(:,:,:,2)),1))) % pres-att vs. pres-unatt
title('ave(T1,T2)')
subplot(3,3,6)
imagesc(squeeze(nanmean((groupMean.PAAU(:,:,:,3)-groupMean.PAAU(:,:,:,4)),1))) % abs-att vs. abs-unatt
% ave(P,A)
subplot(3,3,7)
imagesc(squeeze(nanmean((groupMean.AUT(:,:,:,1,1)-groupMean.AUT(:,:,:,2,1)),1))) % T1-att vs. T1-unatt 
ylabel('ave(P,A)')
subplot(3,3,8)
imagesc(squeeze(nanmean((groupMean.AUT(:,:,:,1,2)-groupMean.AUT(:,:,:,2,2)),1))) % T2-att vs. T2-unatt 
% ave(all)
subplot(3,3,9)
imagesc(squeeze(nanmean((groupMean.AU(:,:,:,1)-groupMean.AU(:,:,:,2)),1))) % att vs. unatt
xlabel('time (s)')
ylabel('frequency (Hz)')
title('ave(all)')
% format subplots
aH = findall(gcf,'type','axes');
paauxtick = [11 61 111];
for iAx = 1:numel(aH)
    axes(aH(iAx));
    rd_timeFreqPlotLabels(twinvals,foi,paauxtick,ytick,0);
    set(gca,'clim',diffClims)
end
colormap(cmap)
rd_supertitle2('amplitude, attended vs. unattended')

% 9 squares, present-absent
fH(5) = figure;
set(gcf,'Position',tf9SquareFigPos)
% T1/T2 x att/unatt
subplot(3,3,1)
imagesc(squeeze(nanmean((groupMean.PAAUT(:,:,:,1,1)-groupMean.PAAUT(:,:,:,3,1)),1))) % T1-pres-att vs. abs-att
ylabel('attended')
title('T1')
subplot(3,3,2)
imagesc(squeeze(nanmean((groupMean.PAAUT(:,:,:,1,2)-groupMean.PAAUT(:,:,:,3,2)),1))) % T2-pres-att vs. abs-att
title('T2')
subplot(3,3,4)
imagesc(squeeze(nanmean((groupMean.PAAUT(:,:,:,2,1)-groupMean.PAAUT(:,:,:,4,1)),1))) % T1-pres-unatt vs. abs-unatt
ylabel('unattended')
subplot(3,3,5)
imagesc(squeeze(nanmean((groupMean.PAAUT(:,:,:,2,2)-groupMean.PAAUT(:,:,:,4,2)),1))) % T2-pres-unatt vs. abs-unatt
% ave(T1,T2)
subplot(3,3,3)
imagesc(squeeze(nanmean((groupMean.PAAU(:,:,:,1)-groupMean.PAAU(:,:,:,3)),1))) % pres-att vs. abs-att
title('ave(T1,T2)')
subplot(3,3,6)
imagesc(squeeze(nanmean((groupMean.PAAU(:,:,:,2)-groupMean.PAAU(:,:,:,4)),1))) % pres-unatt vs. abs-unatt
% ave(A,U)
subplot(3,3,7)
imagesc(squeeze(nanmean((groupMean.PAT(:,:,:,1,1)-groupMean.PAT(:,:,:,2,1)),1))) % T1-pres vs. T1-abs 
ylabel('ave(A,U)')
subplot(3,3,8)
imagesc(squeeze(nanmean((groupMean.PAT(:,:,:,1,2)-groupMean.PAT(:,:,:,2,2)),1))) % T2-pres vs. T2-abs 
% ave(all)
subplot(3,3,9)
imagesc(squeeze(nanmean((groupMean.PA(:,:,:,1)-groupMean.PA(:,:,:,2)),1))) % pres vs. abs
xlabel('time (s)')
ylabel('frequency (Hz)')
title('ave(all)')
% format subplots
aH = findall(gcf,'type','axes');
paauxtick = [11 61 111];
for iAx = 1:numel(aH)
    axes(aH(iAx));
    rd_timeFreqPlotLabels(twinvals,foi,paauxtick,ytick,0);
    set(gca,'clim',diffClims)
end
colormap(cmap)
rd_supertitle2('amplitude, present vs. absent')

%% save figs
if saveFigs
    if strcmp(measure, 'itpc-single-wb')
        fH(fH==0) = [];
        figNames = {'timeFreqSingleByCond','timeFreqSingleAUDiff','timeFreqSinglePADiff'};
    else
        figNames = {'timeFreqSingleByCond','timeFreqSingleAtt','timeFreqSinglePA','timeFreqSingleAUDiff','timeFreqSinglePADiff'};
    end
        
    figPrefix = sprintf('%s_im_wholebrain_%s', figStr, measure);
    rd_saveAllFigs(fH, figNames, figPrefix, figDir)
end

%% Maps
vals = nanmean(nanmean(nanmean(groupMean.amps(:,10:11,toi>0.5 & toi<3.1,1:end-1),4),3),2);
figure
ssm_plotOnMesh(vals', ...
    'alpha, f=[10 11], t=[500 3100]',[], data_hdr, '2d');
colormap(parula)
colorbar
[~,alphaChannels] = sort(vals,'descend');
alphaChannels = alphaChannels(1:5);

vals = nanmean(nanmean(nanmean(groupMean.amps(:,8:14,toi>0.5 & toi<3.1,1:end-1),4),3),2);
figure
ssm_plotOnMesh(vals', ...
    'alpha, f=[8 14], t=[500 3100]',[], data_hdr, '2d');
colormap(parula)
colorbar
[~,alphaChannels] = sort(vals,'descend');
alphaChannels = alphaChannels(1:5);

vals = nanmean(nanmean(nanmean(groupMean.amps(:,40,toi>0.5 & toi<3.1,1:end-1),4),3),2);
figure
ssm_plotOnMesh(vals', ...
    'f=40, t=[500 3100]',[], data_hdr, '2d');
colormap(parula)
colorbar
[~,ssvef40Channels] = sort(vals,'descend');
ssvef40Channels = ssvef40Channels(1:5);

vals = nanmean(nanmean(nanmean(groupMean.amps(:,45:50,toi>0.5 & toi<3.1,1:end-1),4),3),2);
figure
ssm_plotOnMesh(vals', ...
    'gamma, f=[45 50], t=[500 3100]',[], data_hdr, '2d');
colormap(parula)
colorbar
set(gca,'clim',[-.02 .02])

%% select plots
figure
hold on
plot(toi, squeeze(mean(mean(mean(groupMean.amps(alphaChannels,10:11,:,1:end-1),4),2),1)))
plot(toi, squeeze(mean(mean(mean(groupMean.amps(ssvef40Channels,40,:,1:end-1),4),2),1)), 'r')

figure
plot(toi, squeeze(mean(mean(groupMean.amps(alphaChannels,10:11,:,:),2),1)))

figure
plot(toi, squeeze(mean(mean(groupMean.ampsAtt(alphaChannels,10:11,:,:),2),1)))
legend('att T1','att T2')

figure
plot(toi, squeeze(mean(mean(groupMean.ampsPA(alphaChannels,10:11,:,:),2),1)))
legend('PP','AP','PA','AA')

figure
plot(toi, squeeze(mean(mean(groupMean.ampsPA(alphaChannels,8:14,:,:),2),1)))
legend('PP','AP','PA','AA')

%% multiplots
% channel plot setup
cfg = [];
cfg.layout = layout;
cfg.colormap = cmap;
cfg.zlim = [-5 5];

TFdata.label = data_hdr.label(1:numel(A.channels));
TFdata.dimord = 'chan_freq_time';
TFdata.freq = foi;
TFdata.time = twinvals;

% A-U amp
% vals = groupMean.AU(:,:,:,1)-groupMean.AU(:,:,:,2); % chan x freq x time
% TFdata.powspctrm = vals;
% cfg.zlim = diffClims;
% 
% fH = [];
% figPos = [1 1 1250 930];
% fH(1) = figure('Position',figPos);
% ft_multiplotTFR(cfg, TFdata);
% rd_supertitle2('amplitude, A-U')

% A-U stats
vals = groupTStat.AU; % chan x freq x time
TFdata.powspctrm = vals;

fH = [];
figPos = [1 1 1250 930];
fH(1) = figure('Position',figPos);
ft_multiplotTFR(cfg, TFdata);
rd_supertitle2('A-U t-stat')

% P-A stats
vals = groupTStat.PA; % chan x freq x time
TFdata.powspctrm = vals;

figPos = [1 1 1250 930];
fH(2) = figure('Position',figPos);
ft_multiplotTFR(cfg, TFdata);
rd_supertitle2('P-A t-stat')

% A-U stats T1
for iT = 1:2
    vals = groupTStat.AUT(:,:,:,iT); % chan x freq x time
    TFdata.powspctrm = vals;
    
    figPos = [1 1 1250 930];
    fH(3+iT-1) = figure('Position',figPos);
    ft_multiplotTFR(cfg, TFdata);
    rd_supertitle2(sprintf('A-U T%d t-stat',iT))
end

if saveFigs
    figPrefix = sprintf('%s_immap_wholebrain_%s', figStr, measure);
    rd_saveAllFigs(fH, {'timeFreqSingleAUDiffTStat','timeFreqSinglePADiffTStat','timeFreqSingleAUDiffT1TStat','timeFreqSingleAUDiffT2TStat'}, figPrefix, figDir)
end

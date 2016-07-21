function rd_plotTADetectDiscrimGroupTimeFreqSingle(A, measure, groupMean, saveFigs, figDir, figStr) 

%% args
if nargin<4
    saveFigs = 0;
end
if nargin<6
    figStr = '';
    if saveFigs==1
        error('If youre saving figs, you should specify a figStr')
    end
end

figTitle = und2space(figStr);

%% setup
plotOrder = [1 5 3 7 2 6 4 8 9];

tf9SquareFigPos = [50 50 850 850];
tf6SquareFigPos = [50 50 850 530];

cmap = flipud(lbmap(64,'RedBlue'));
% cmap = colormap;

foi = A.stfFoi;
twinvals = A.stftwinvals;

ytick = 10:10:numel(foi);
paauxtick = [11 61 111];

switch measure
    case 'stf-single'
        clims = [-0.1 0.1];
        diffClims = [-0.07 0.07];
    otherwise
        error('measure not recognized')
end

%% 9 squares, attended-unattended
fH = [];
fH(1) = figure;
set(gcf,'Position',tf9SquareFigPos)
% T1/T2 x pres/abs
subplot(3,3,1)
imagesc(groupMean.PAAUT(:,:,1,1)-groupMean.PAAUT(:,:,2,1)) % T1-pres-att vs. unatt
ylabel('present')
title('T1')
subplot(3,3,2)
imagesc(groupMean.PAAUT(:,:,1,2)-groupMean.PAAUT(:,:,2,2)) % T2-pres-att vs. unatt
title('T2')
subplot(3,3,4)
imagesc(groupMean.PAAUT(:,:,3,1)-groupMean.PAAUT(:,:,4,1)) % T1-abs-att vs. unatt
ylabel('absent')
subplot(3,3,5)
imagesc(groupMean.PAAUT(:,:,3,2)-groupMean.PAAUT(:,:,4,2)) % T2-abs-att vs. unatt
% ave(T1,T2)
subplot(3,3,3)
imagesc(groupMean.PAAU(:,:,1)-groupMean.PAAU(:,:,2)) % pres-att vs. pres-unatt
title('ave(T1,T2)')
subplot(3,3,6)
imagesc(groupMean.PAAU(:,:,3)-groupMean.PAAU(:,:,4)) % abs-att vs. abs-unatt
% ave(P,A)
subplot(3,3,7)
imagesc(groupMean.AUT(:,:,1,1)-groupMean.AUT(:,:,2,1)) % T1-att vs. T1-unatt 
ylabel('ave(P,A)')
subplot(3,3,8)
imagesc(groupMean.AUT(:,:,1,2)-groupMean.AUT(:,:,2,2)) % T2-att vs. T2-unatt 
% ave(all)
subplot(3,3,9)
imagesc(groupMean.AU(:,:,1)-groupMean.AU(:,:,2)) % att vs. unatt
xlabel('time (s)')
ylabel('frequency (Hz)')
title('ave(all)')
% format subplots
aH = findall(gcf,'type','axes');
for iAx = 1:numel(aH)
    axes(aH(iAx));
    rd_timeFreqPlotLabels(twinvals,foi,paauxtick,ytick,0);
    set(gca,'clim',diffClims)
end
colormap(cmap)
rd_supertitle2('attended vs. unattended')

%% 9 squares, present-absent
fH(2) = figure;
set(gcf,'Position',tf9SquareFigPos)
% T1/T2 x att/unatt
subplot(3,3,1)
imagesc(groupMean.PAAUT(:,:,1,1)-groupMean.PAAUT(:,:,3,1)) % T1-pres-att vs. abs-att
ylabel('attended')
title('T1')
subplot(3,3,2)
imagesc(groupMean.PAAUT(:,:,1,2)-groupMean.PAAUT(:,:,3,2)) % T2-pres-att vs. abs-att
title('T2')
subplot(3,3,4)
imagesc(groupMean.PAAUT(:,:,2,1)-groupMean.PAAUT(:,:,4,1)) % T1-pres-unatt vs. abs-unatt
ylabel('unattended')
subplot(3,3,5)
imagesc(groupMean.PAAUT(:,:,2,2)-groupMean.PAAUT(:,:,4,2)) % T2-pres-unatt vs. abs-unatt
% ave(T1,T2)
subplot(3,3,3)
imagesc(groupMean.PAAU(:,:,1)-groupMean.PAAU(:,:,3)) % pres-att vs. abs-att
title('ave(T1,T2)')
subplot(3,3,6)
imagesc(groupMean.PAAU(:,:,2)-groupMean.PAAU(:,:,4)) % pres-unatt vs. abs-unatt
% ave(A,U)
subplot(3,3,7)
imagesc(groupMean.PAT(:,:,1,1)-groupMean.PAT(:,:,2,1)) % T1-pres vs. T1-abs 
ylabel('ave(A,U)')
subplot(3,3,8)
imagesc(groupMean.PAT(:,:,1,2)-groupMean.PAT(:,:,2,2)) % T2-pres vs. T2-abs 
% ave(all)
subplot(3,3,9)
imagesc(groupMean.PA(:,:,1)-groupMean.PA(:,:,2)) % pres vs. abs
xlabel('time (s)')
ylabel('frequency (Hz)')
title('ave(all)')
% format subplots
aH = findall(gcf,'type','axes');
for iAx = 1:numel(aH)
    axes(aH(iAx));
    rd_timeFreqPlotLabels(twinvals,foi,paauxtick,ytick,0);
    set(gca,'clim',diffClims)
end
colormap(cmap)
rd_supertitle2('present vs. absent')

%% 6 squares, present
fH(3) = figure;
set(gcf,'Position',tf6SquareFigPos)
% T1/T2 x att/unatt
subplot(2,3,1)
imagesc(groupMean.PAAUT(:,:,1,1)) % T1-pres-att
ylabel('attended')
title('T1')
subplot(2,3,2)
imagesc(groupMean.PAAUT(:,:,1,2)) % T2-pres-att
title('T2')
subplot(2,3,4)
imagesc(groupMean.PAAUT(:,:,2,1)) % T1-pres-unatt
ylabel('unattended')
subplot(2,3,5)
imagesc(groupMean.PAAUT(:,:,2,2)) % T2-pres-unatt
% ave(T1,T2)
subplot(2,3,3)
imagesc(groupMean.PAAU(:,:,1)) % pres-att
title('ave(T1,T2)')
subplot(2,3,6)
imagesc(groupMean.PAAU(:,:,2)) % pres-unatt
xlabel('time (s)')
ylabel('frequency (Hz)')
title('ave(all)')
% format subplots
aH = findall(gcf,'type','axes');
for iAx = 1:numel(aH)
    axes(aH(iAx));
    rd_timeFreqPlotLabels(twinvals,foi,paauxtick,ytick,0);
    set(gca,'clim',clims)
end
colormap(cmap)
rd_supertitle2('target present')

%% 6 squares, absent
fH(4) = figure;
set(gcf,'Position',tf6SquareFigPos)
% T1/T2 x att/unatt
subplot(2,3,1)
imagesc(groupMean.PAAUT(:,:,3,1)) % T1-abs-att
ylabel('attended')
title('T1')
subplot(2,3,2)
imagesc(groupMean.PAAUT(:,:,3,2)) % T2-abs-att
title('T2')
subplot(2,3,4)
imagesc(groupMean.PAAUT(:,:,4,1)) % T1-abs-unatt
ylabel('unattended')
subplot(2,3,5)
imagesc(groupMean.PAAUT(:,:,4,2)) % T2-abs-unatt
% ave(T1,T2)
subplot(2,3,3)
imagesc(groupMean.PAAU(:,:,3)) % abs-att
title('ave(T1,T2)')
subplot(2,3,6)
imagesc(groupMean.PAAU(:,:,4)) % abs-unatt
xlabel('time (s)')
ylabel('frequency (Hz)')
title('ave(all)')
% format subplots
aH = findall(gcf,'type','axes');
for iAx = 1:numel(aH)
    axes(aH(iAx));
    rd_timeFreqPlotLabels(twinvals,foi,paauxtick,ytick,0);
    set(gca,'clim',clims)
end
colormap(cmap)
rd_supertitle2('target absent')

%% save figs
if saveFigs
    figPrefix = sprintf('%s_im', figStr);
    switch measure
        case 'stf-single'
            figNames = {'timeFreqSingleAttVsUnatt','timeFreqSinglePVsA','timeFreqSinglePres','timeFreqSingleAbs'};
        otherwise
            error('measure not recognized')
    end
    rd_saveAllFigs(fH, figNames, figPrefix, figDir)
end
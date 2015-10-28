function rd_plotTADetectDiscrimGroupTimeFreq(A, measure, groupMean, saveFigs, figDir, figStr) 

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

tf9FigPos = [0 250 1280 580];
tf3FigPos = [200 475 1000 275];

foi = A.tfFoi;
toi = A.tfToi;
eventTimes = A.eventTimes;
trigNames = A.trigNames;
attNames = A.attNames;
paNames = A.PANames;

tfAmps = groupMean.amps;
tfAmpsAtt = groupMean.ampsAtt;
tfAmpsPA = groupMean.ampsPA;
tfPADiff = groupMean.paDiff;

t1PADiff = tfPADiff(:,:,1);
t2PADiff = tfPADiff(:,:,2);
nTrigs = numel(trigNames);

ytick = 10:10:numel(foi);
xtick = 51:50:numel(toi);
hack = plotOrder;
hack(hack>4) = hack(hack>4)+1;

switch measure
    case 'tf'
        clims = [0 18];
        diffClims = [-4 4];
    case 'stf'
        clims = [0 50];
        diffClims = [-5 5];
    otherwise
        error('measure not recognized')
end

%% fig 1
fH = [];
fH(1) = figure;
set(gcf,'Position',tf9FigPos)
for iTrig = 1:nTrigs
    subplot(2,5,hack(iTrig))
    imagesc(tfAmps(:,:,iTrig),clims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    if iTrig==nTrigs
        xlabel('time (s)')
        ylabel('frequency (Hz)')
    end
    title(trigNames{iTrig})
end
rd_supertitle(figTitle);
rd_raiseAxis(gca);

%% fig 2
fH(2) = figure;
set(gcf,'Position',tf3FigPos)
for iAtt = 1:size(tfAmpsAtt,3)
    subplot(1,3,iAtt)
    imagesc(tfAmpsAtt(:,:,iAtt),clims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    xlabel('time (s)')
    ylabel('frequency (Hz)')
    title(attNames{iAtt})
end
subplot(1,3,3)
imagesc(tfAmpsAtt(:,:,2)-tfAmpsAtt(:,:,1),diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (s)')
ylabel('frequency (Hz)')
title('attT2 - attT1')
rd_supertitle(figTitle);
rd_raiseAxis(gca);

%% fig 3
fH(3) = figure;
set(gcf,'Position',tf9FigPos)
for iPA = 1:size(tfAmpsPA,3)
    subplot(2,4,iPA)
    imagesc(tfAmpsPA(:,:,iPA),clims)
    rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
    xlabel('time (s)')
    ylabel('frequency (Hz)')
    title(paNames{iPA})
end
subplot(2,4,5)
imagesc(t1PADiff,diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (s)')
ylabel('frequency (Hz)')
title('T1 P-A')
subplot(2,4,6)
imagesc(t2PADiff,diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (s)')
ylabel('frequency (Hz)')
title('T2 P-A')
subplot(2,4,7)
imagesc(t2PADiff - t1PADiff,diffClims)
rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
xlabel('time (s)')
ylabel('frequency (Hz)')
title('T2 vs. T1 P-A')
rd_supertitle(figTitle);
rd_raiseAxis(gca);

if saveFigs
    figPrefix = sprintf('%s_im', figStr);
    switch measure
        case 'tf'
            figNames = {'timeFreqByCond','timeFreqAtt','timeFreqPA'};
        case 'stf'
            figNames = {'timeFreqSingleByCond','timeFreqSingleAtt','timeFreqSinglePA'};
        otherwise
            error('measure not recognized')
    end
    rd_saveAllFigs(fH, figNames, figPrefix, figDir)
end

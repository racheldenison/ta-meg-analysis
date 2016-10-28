function rd_plotTADetectDiscrimGroupAmpsWholebrain(A, measure, subjects, ...
    groupData, groupMean, groupSte, ...
    saveFigs, figDir, figStr)

%% setup
load data/data_hdr.mat
twin = A.wtwin;
wPAAUT = groupMean.PAAUT;
wPAAU = groupMean.PAAU;
wAUT = groupMean.AUT;
wPAT = groupMean.PAT;

%% plot setup
paauNames = {'P-att','P-unatt','A-att','A-unatt'};
auNames = {'att','unatt'};
paNames = {'P','A'};
switch A.normalizeOption
    case 'none'
        clims = [0 200];
        diffClims = [-40 40];
    case 'stim'
        clims = [0.5 1.5];
        diffClims = [-.3 .3];
end
cmap = flipud(lbmap(64,'RedBlue'));
twindow = twin(1):twin(end);
nBins = 6;
binSize = round(numel(twindow)/nBins);
load parula

%% movie
tstep = 10;
iPAAU = 1;
iT = 1;
figure
for iTime = 1:tstep:numel(twindow)
    ssm_plotOnMesh(wPAAUT(iTime,:,iPAAU,iT), ...
        sprintf('wPAAUT t=%d',twindow(iTime)),[], data_hdr, '2d');
    set(gca,'CLim',clims)
    colormap(parula)
    pause(.1)
end

%% time bins
fH = [];
% PAAU
figPos = [32 150 200*nBins 750];
for iT = 1:2
    fH(1+iT-1) = figure('Position',figPos);
    for iPAAU = 1:4
        for iBin = 1:nBins
            subplot(4,nBins,iBin + nBins*(iPAAU-1))
            tidx = (1:binSize+1) + (iBin-1)*binSize;
            str = sprintf('wPAAUT, %s, t=[%d %d]',paauNames{iPAAU}, twindow(tidx(1)), twindow(tidx(end)));
            vals = mean(wPAAUT(tidx,:,iPAAU,iT),1);
            ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
            set(gca,'CLim',clims)
            colormap(parula)
        end
    end
    rd_supertitle2(sprintf('T%d',iT))
end

% AU
figPos = [32 250 200*nBins 650];
for iT = 1:2
    fH(3+iT-1) = figure('Position',figPos);
    for iAU = 1:2
        for iBin = 1:nBins
            subplot(3,nBins,iBin + nBins*(iAU-1))
            tidx = (1:binSize+1) + (iBin-1)*binSize;
            str = sprintf('wAUT, %s, t=[%d %d]',auNames{iAU}, twindow(tidx(1)), twindow(tidx(end)));
            vals = mean(wAUT(tidx,:,iAU,iT),1);
            ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
            set(gca,'CLim',clims)
            colormap(parula)
            freezeColors
        end
    end
    for iBin = 1:nBins
        subplot(3,nBins,iBin + nBins*2)
        tidx = (1:binSize+1) + (iBin-1)*binSize;
        str = sprintf('wAUT, A-U, t=[%d %d]',twindow(tidx(1)), twindow(tidx(end)));
        vals = mean((wAUT(tidx,:,1,iT) - wAUT(tidx,:,2,iT)),1);
        ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
        set(gca,'CLim',diffClims)
    end
    colormap(cmap)
    rd_supertitle2(sprintf('T%d',iT))
end

% PA
figPos = [32 250 200*nBins 650];
for iT = 1:2
    fH(5+iT-1) = figure('Position',figPos);
    for iPA = 1:2
        for iBin = 1:nBins
            subplot(3,nBins,iBin + nBins*(iPA-1))
            tidx = (1:binSize+1) + (iBin-1)*binSize;
            str = sprintf('wPAT, %s, t=[%d %d]',paNames{iPA}, twindow(tidx(1)), twindow(tidx(end)));
            vals = mean(wPAT(tidx,:,iPA,iT),1);
            ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
            set(gca,'CLim',clims)
            colormap(parula)
            freezeColors
        end
    end
    for iBin = 1:nBins
        subplot(3,nBins,iBin + nBins*2)
        tidx = (1:binSize+1) + (iBin-1)*binSize;
        str = sprintf('wPAT, P-A, t=[%d %d]',twindow(tidx(1)), twindow(tidx(end)));
        vals = mean((wPAT(tidx,:,1,iT) - wPAT(tidx,:,2,iT)),1);
        ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
        set(gca,'CLim',diffClims)
    end
    colormap(cmap)
    rd_supertitle2(sprintf('T%d',iT))
end

% AUDiff for present and absent separately
figPos = [32 250 200*nBins 450];
for iT = 1:2
    fH(7+iT-1) = figure('Position',figPos);
    for iBin = 1:nBins
        subplot(2,nBins,iBin)
        tidx = (1:binSize+1) + (iBin-1)*binSize;
        str = sprintf('%s, t=[%d %d]','P-att - P-unatt', twindow(tidx(1)), twindow(tidx(end)));
        vals = mean(wPAAUT(tidx,:,1,iT),1) - mean(wPAAUT(tidx,:,2,iT),1);
        ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
        set(gca,'CLim',diffClims)
        colormap(cmap)
    end
    for iBin = 1:nBins
        subplot(2,nBins,iBin + nBins)
        tidx = (1:binSize+1) + (iBin-1)*binSize;
        str = sprintf('%s, t=[%d %d]','A-att - A-unatt', twindow(tidx(1)), twindow(tidx(end)));
        vals = mean(wPAAUT(tidx,:,3,iT),1) - mean(wPAAUT(tidx,:,4,iT),1);
        ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
        set(gca,'CLim',diffClims)
        colormap(cmap)
    end
    rd_supertitle2(sprintf('T%d',iT))
end

% AUDiff for present and absent separately, T1 & T2 combined
figPos = [32 250 200*nBins 450];
fH(9) = figure('Position',figPos);
for iBin = 1:nBins
    subplot(2,nBins,iBin)
    tidx = (1:binSize+1) + (iBin-1)*binSize;
    str = sprintf('%s, t=[%d %d]','P-att - P-unatt', twindow(tidx(1)), twindow(tidx(end)));
    vals = mean(wPAAU(tidx,:,1),1) - mean(wPAAU(tidx,:,2),1);
    ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
    set(gca,'CLim',diffClims)
    colormap(cmap)
end
for iBin = 1:nBins
    subplot(2,nBins,iBin + nBins)
    tidx = (1:binSize+1) + (iBin-1)*binSize;
    str = sprintf('%s, t=[%d %d]','A-att - A-unatt', twindow(tidx(1)), twindow(tidx(end)));
    vals = mean(wPAAU(tidx,:,3),1) - mean(wPAAU(tidx,:,4),1);
    ssm_plotOnMesh(vals,str,[], data_hdr, '2d');
    set(gca,'CLim',diffClims)
    colormap(cmap)
end
rd_supertitle2('T1 & T2')

if saveFigs
    figPrefix = sprintf('%s_map_wholebrain_%dHz', figStr, ssvefFreq);
    rd_saveAllFigs(fH, {'wPAAUT1','wPAAUT2','wAUT1','wAUT2','wPAT1','wPAT2','wPAAUDiffT1','wPAAUDiffT2','wPAAUDiffT1T2Comb'}, figPrefix, figDir)
end


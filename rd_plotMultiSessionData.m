% rd_plotMultiSessionData.m

%% Setup
exptType = 'TANoise';
exptDir = '/Local/Users/denison/Data/TANoise/MEG';
% sessionDirs = {'R0817_20171212','R0817_20171213'};
% sessionDirs = {'R1187_20180105','R1187_20180108'};
sessionDirs = {'R0983_20180111','R0983_20180112'};
trialsOption = 'trialAve'; % 'singleTrials','trialAve'
wIdx = 1:6701;
tfIdx = 1:671;

analStr = 'ebi'; % '', 'ebi', etc.
ssvefFreq = 20;
nTopChannels = 5; % 1, 5, etc., or [] for iqrThresh
% iqrThresh = []; % 10, or [] for nTopChannels
% weightChannels = 0; % weight channels according to average SSVEF amp - only works for top channels
trialSelection = 'all'; % 'all','validCorrect', etc
respTargetSelection = ''; % '','T1Resp','T2Resp'
excludeTrialsFt = 1;

nSessions = numel(sessionDirs);
channelSelectionStr = sprintf('topChannels%d', nTopChannels);
if excludeTrialsFt
    analStr = sprintf('%s_ft', analStr);
end
switch trialsOption
    case 'trialAve'
        trialsStr = '';
    case 'singleTrials'
        trialsStr = '_singleTrials'; 
    otherwise
        error('trialsOption not recognized')
end

%% Load data
for iSession = 1:nSessions
    sessionDir = sessionDirs{iSession};
    dataDir = sprintf('%s/%s', exptDir, sessionDir);
    matDir = sprintf('%s/mat', dataDir);
    
    switch analStr
        case ''
            dataFile = dir(sprintf('%s/analysis%s_R*_%s_%sTrials%s_%dHz.mat', matDir, trialsStr, channelSelectionStr, trialSelection, respTargetSelection, ssvefFreq));
        otherwise
            dataFile = dir(sprintf('%s/analysis%s_R*_%s_%s_%sTrials%s_%dHz.mat', matDir, trialsStr, analStr, channelSelectionStr, trialSelection, respTargetSelection, ssvefFreq));
    end
    
    a0 = load(sprintf('%s/%s', matDir, dataFile.name));
    A(iSession) = a0.A;
end

nA = numel(A);

trigNames = A(1).trigNames;
nTrigs = numel(trigNames);
t = A(2).t;
eventTimes = A(1).eventTimes;
attNames = A(1).attNames;

%% Plotting setup
plotOrder = [1 5 3 7 2 6 4 8 9];
extendedMap = flipud(lbmap(nTrigs-1+4,'RedBlue'));
selectedMap = extendedMap([1:(nTrigs-1)/2 (end-(nTrigs-1)/2)+1:end],:);
trigColors = [selectedMap; 0 0 0];
trigBlue = mean(selectedMap(1:(nTrigs-1)/2,:));
trigRed = mean(selectedMap((end-(nTrigs-1)/2)+1:end,:));

tsFigPos = [0 500 1250 375];
% ts2FigPos = [0 500 1100 600];
% ts3FigPos = [0 500 1100 900];
% condFigPos = [250 300 750 650];
% tf9FigPos = [0 250 1280 580];
tf3FigPos = [200 475 980 330];

set(0,'defaultLineLineWidth',1)

switch exptType
    case 'TADetectDiscrim'
        xtickint = 50;
    case 'TAContrast'
        xtickint = 100;
    case 'TANoise'
        xtickint = 100;
end

%% Combine and plot
% different analyses for trialAve and singleTrials
switch trialsOption
    case 'trialAve'
        %% wAmps
        vals = [];
        for iA = 1:nA
            vals(:,:,iA) = A(iA).wAmps(wIdx,:);
        end
        wAmps = nanmean(vals, 3);
        
        fH = [];
        fH(1) = figure;
        set(gcf,'Position',tsFigPos)
        set(gca,'ColorOrder',trigColors)
        hold all
        % plot(t, wAmps(:,plotOrder))
        plot(t, wAmps(:,end), 'k') % blank
        plot(t, nanmean(wAmps(:,plotOrder(1:(nTrigs-1)/2)),2),'color',trigBlue,'LineWidth',4)
        plot(t, nanmean(wAmps(:,plotOrder(end-(nTrigs-1)/2):end-1),2),'color',trigRed,'LineWidth',4)
        for iEv = 1:numel(eventTimes)
            vline(eventTimes(iEv),'k');
        end
        % legend(trigNames(plotOrder))
        legend('blank','att T1','att T2')
        xlabel('time (ms)')
        ylabel('wavelet amp')
        % title([sprintf('%d Hz, channel', ssvefFreq) sprintf(' %d', channels) wstrt])
        
        %% time freq single
        vals = [];
        for iA = 1:nA
            vals(:,:,:,iA) = A(iA).stfAmpsAtt(:,tfIdx,:);
        end
        tfSingleAmpsAtt = nanmean(vals, 4);
        tfSingleAmpsAttDiff = tfSingleAmpsAtt(:,:,2)-tfSingleAmpsAtt(:,:,1);
        maxval = max(tfSingleAmpsAtt(:));
        maxvaldiff = max(abs(tfSingleAmpsAttDiff(:)));
        
        % figures
        toi = A(1).stfToi;
        foi = A(2).stfFoi;
        ytick = 10:10:numel(foi);
        xtick = 51:xtickint:numel(toi);
        clims = [0 maxval];
        diffClims = [-maxvaldiff maxvaldiff];
        
        fH(2) = figure;
        set(gcf,'Position',tf3FigPos)
        attNames = {'attT1','attT2'};
        for iAtt = 1:size(tfSingleAmpsAtt,3)
            subplot(1,3,iAtt)
            imagesc(tfSingleAmpsAtt(:,:,iAtt))
            imagesc(tfSingleAmpsAtt(:,:,iAtt),clims)
            rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
            xlabel('time (s)')
            ylabel('frequency (Hz)')
            title(attNames{iAtt})
        end
        subplot(1,3,3)
        imagesc(tfSingleAmpsAttDiff,diffClims)
        rd_timeFreqPlotLabels(toi,foi,xtick,ytick,eventTimes);
        xlabel('time (s)')
        ylabel('frequency (Hz)')
        title('attT2 - attT1')
        % rd_supertitle(['channel' sprintf(' %d', channels) wstrt]);
        rd_raiseAxis(gca);
        
    case 'singleTrials'
        %% wAmps single
        wAmpsAtt = cat(2, A(1).wAmpsAtt(wIdx,:,:), A(2).wAmpsAtt(wIdx,:,:));
        
        fH(3) = figure;
        set(gcf,'Position',tsFigPos)
        hold on
        plot(t, nanmean(wAmpsAtt(:,:,1),2),'color',trigBlue,'LineWidth',4)
        plot(t, nanmean(wAmpsAtt(:,:,2),2),'color',trigRed,'LineWidth',4)
        legend(attNames)
        [~, emp, err] = rd_bootstrapCI(wAmpsAtt(:,:,1)');
        shadedErrorBar(t, emp, err, {'color',trigBlue,'LineWidth',4}, 1)
        [~, emp, err] = rd_bootstrapCI(wAmpsAtt(:,:,2)');
        shadedErrorBar(t, emp, err, {'color',trigRed,'LineWidth',4}, 1)
        for iEv = 1:numel(eventTimes)
            vline(eventTimes(iEv),'k');
        end
        xlabel('time (ms)')
        ylabel('single trial wavelet amp')
        
        %% itpc 
        % Method 1: average sessions without recomputing ITPC
%         vals = [];
%         for iA = 1:nA
%             vals(:,:,iA) = A(iA).wITPCAtt;
%         end
%         wITPCAtt = nanmean(vals, 3);

        % Method 2: recompute ITPC across all trials
%         % combine spectra for all trials
%         vals = [];
%         for iA = 1:nA
%             vals = cat(2, vals, A(iA).wSpecAtt);
%         end
% 
%         % define itpc function
%         itpcFun = @(spectrum) squeeze(abs(nanmean(exp(1i*angle(spectrum)),1))); % mean across trials
%         % generate sample indices, same for all conditions
%         nBoot = 100;
%         [~, bootIdx] = bootstrp(nBoot, @(x) [], vals(:,:,1,1)');
%         % recompute itpc from all trials in condition
%         fprintf('\nbootstrap itpc:\n')
%         for iAtt = 1:numel(attNames)
%             fprintf('%s\n', attNames{iAtt})
%             for iCh = 1:nTopChannels
%                 fprintf('ch %d\n', iCh)
%                 spectrum = vals(:,:,iAtt,iCh)';
%                 emp = itpcFun(spectrum);
%                 
%                 itpcResamp = [];
%                 for iBoot = 1:nBoot
%                     itpcResamp(:,iBoot) = itpcFun(spectrum(bootIdx(:,iBoot),:));
%                 end
% 
%                 wITPCAtt0(:,iCh,iAtt) = emp;
%                 wITPCAttResamp0(:,:,iCh,iAtt) = itpcResamp;
%             end
%         end
%         
%         % mean across channels
%         chIdx = 1:nTopChannels;
%         wITPCAtt = squeeze(nanmean(wITPCAtt0(:,chIdx,:),2)); 
%         wITPCAttResamp = squeeze(nanmean(wITPCAttResamp0(:,:,chIdx,:),3));
%         
%         % ci and error bars
%         wITPCAttCI = prctile(wITPCAttResamp, [2.5 97.5], 2);
%         for iAtt = 1:numel(attNames)
%             wITPCAttErr(:,1,iAtt) = wITPCAttCI(:,1,iAtt)-wITPCAtt(:,iAtt);
%             wITPCAttErr(:,2,iAtt) = wITPCAtt(:,iAtt)-wITPCAttCI(:,2,iAtt);
%         end
        
        % Method 3: recompute separately for each session
        % define itpc function
        itpcFun = @(spectrum) squeeze(abs(nanmean(exp(1i*angle(spectrum)),1))); % mean across trials
        nBoot = 100;
        
        % combine spectra for all trials
        fprintf('\nbootstrap itpc\n')
        for iA = 1:nA
            vals = A(iA).wSpecAtt(wIdx,:,:,:);
            % generate sample indices, same for all conditions
            [~, bootIdx] = bootstrp(nBoot, @(x) [], vals(:,:,1,1)');
            % recompute itpc from all trials in condition
            fprintf('\nsession %d\n', iA)
            for iAtt = 1:numel(attNames)
                fprintf('%s\n', attNames{iAtt})
                for iCh = 1:nTopChannels
                    fprintf('ch %d\n', iCh)
                    spectrum = vals(:,:,iAtt,iCh)';
                    emp = itpcFun(spectrum);
                    
                    itpcResamp = [];
                    for iBoot = 1:nBoot
                        itpcResamp(:,iBoot) = itpcFun(spectrum(bootIdx(:,iBoot),:));
                    end
                    
                    wITPCAtt0(:,iCh,iAtt,iA) = emp;
                    wITPCAttResamp0(:,:,iCh,iAtt,iA) = itpcResamp;
                end
            end
        end
        
        % mean across channels and sessions
        chIdx = 1:nTopChannels;
        wITPCAtt = squeeze(nanmean(nanmean(wITPCAtt0(:,chIdx,:,:),2),4)); 
        wITPCAttResamp = squeeze(nanmean(nanmean(wITPCAttResamp0(:,:,chIdx,:,:),3),5));
        
        % ci and error bars
        wITPCAttCI = prctile(wITPCAttResamp, [2.5 97.5], 2);
        for iAtt = 1:numel(attNames)
            wITPCAttErr(:,1,iAtt) = wITPCAttCI(:,1,iAtt)-wITPCAtt(:,iAtt);
            wITPCAttErr(:,2,iAtt) = wITPCAtt(:,iAtt)-wITPCAttCI(:,2,iAtt);
        end
        

        fH(4) = figure;
        set(gcf,'Position',tsFigPos)
        hold on
        plot(t, wITPCAtt(:,1),'color',trigBlue,'LineWidth',4)
        plot(t, wITPCAtt(:,2),'color',trigRed,'LineWidth',4)
        legend(attNames)
        shadedErrorBar(t, wITPCAtt(:,1), wITPCAttErr(:,:,1), {'color',trigBlue,'LineWidth',4}, 1)
        shadedErrorBar(t, wITPCAtt(:,2), wITPCAttErr(:,:,2), {'color',trigRed,'LineWidth',4}, 1)
        legend(attNames)
        for iEv = 1:numel(eventTimes)
            vline(eventTimes(iEv),'k');
        end
        xlabel('time (ms)')
        ylabel('wavelet itpc')    
end



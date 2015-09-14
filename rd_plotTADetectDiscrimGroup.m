% rd_plotTADetectDiscrimGroup.m

%% Setup
exptDir = '/Volumes/DRIVE1/DATA/rachel/MEG/TADetectDiscrim/MEG';
analStr = 'ebi_ft'; % '', 'ebi', etc.
ssvefFreq = 30;
nTopChannels = 5; % 1, 5, etc.
measure = 'w'; % wAmpsAtt wAmpsPA hAmpsAtt hAmpsPA tf stf w

subjects = {'R0817_20150504', 'R0973_20150727', 'R0974_20150728', ...
    'R0861_20150813', ...
    'R0504_20150805', 'R0983_20150813', 'R0898_20150828','R0436_20150904'};
% subjects = {'R0817_20150504', 'R0973_20150727', ...
%     'R0861_20150813', 'R0898_20150828'};
% subjects = {'R0973_20150727', 'R0974_20150728', ...
%     'R0861_20150813', ...
%     'R0504_20150805', 'R0983_20150813', 'R0898_20150828','R0436_20150904'};

nSubjects = numel(subjects);

tstart = -500; % ms
tstop = 3600; % ms
t = tstart:tstop;
evTimes = [0 500 1500 2100 3100];
eventTimes = evTimes;

%% Get data
for iSubject = 1:nSubjects
    subject = subjects{iSubject};
    
    sessionDir = subject;
    dataDir = sprintf('%s/%s', exptDir, sessionDir);
    matDir = sprintf('%s/mat', dataDir);
    
    analysisFile = dir(sprintf('%s/analysis_*_%s_topChannels%d_%dHz.mat', matDir, analStr, nTopChannels, ssvefFreq));
    
    if numel(analysisFile)==1
        load(sprintf('%s/%s', matDir, analysisFile.name))
    else
        error('too many or too few matching analysis files')
    end
    
    switch measure
        case 'w'
            groupData.amps(:,:,iSubject) = A.wAmps;
            groupData.ampsAtt(:,:,iSubject) = A.wAmpsAtt;
            groupData.ampsPA(:,:,iSubject) = A.wAmpsPA;
        case 'h'
            groupData.amps(:,:,iSubject) = A.hAmps;
            groupData.ampsAtt(:,:,iSubject) = A.hAmpsAtt;
            groupData.ampsPA(:,:,iSubject) = A.hAmpsPA;
        case 'tf'
            groupData.amps(:,:,:,iSubject) = A.tfAmps;
            groupData.ampsAtt(:,:,:,iSubject) = A.tfAmpsAtt;
            groupData.ampsPA(:,:,:,iSubject) = A.tfAmpsPA;
            groupData.paDiff(:,:,:,iSubject) = A.tfPADiff;
        case 'stf'
            groupData.amps(:,:,:,iSubject) = A.stfAmps;
            groupData.ampsAtt(:,:,:,iSubject) = A.stfAmpsAtt;
            groupData.ampsPA(:,:,:,iSubject) = A.stfAmpsPA;
            groupData.paDiff(:,:,:,iSubject) = A.stfPADiff;
        otherwise
            v = A.(measure);
            
            if numel(size(v))==2
                vals(:,:,iSubject) = v;
            elseif numel(size(v))==3
                vals(:,:,:,iSubject) = v;
            else
                error('too many or too few dimensions')
            end
    end
end

if exist('vals','var')
    dataType = 'vals';
else
    dataType = 'groupData';
end

%% Plots
switch dataType
    case 'vals'
        %% Calculate
        valsDiff = diff(vals);
        
        %% summary stats
        sdim = numel(size(vals)); % subject dimension
        
        valsMean = squeeze(mean(vals,sdim));
        valsSte = squeeze(std(vals,0,sdim)./sqrt(nSubjects));
        
        valsDiffMean = squeeze(mean(valsDiff,sdim));
        valsDiffSte = squeeze(std(valsDiff,0,sdim)./sqrt(nSubjects));
        
        %% plot setup
        nTrigs = numel(A.trigNames);
        plotOrder = [1 5 3 7 2 6 4 8 9];
        extendedMap = flipud(lbmap(nTrigs-1+4,'RedBlue'));
        selectedMap = extendedMap([1:(nTrigs-1)/2 (end-(nTrigs-1)/2)+1:end],:);
        trigColors = [selectedMap; 0 0 0];
        trigBlue = mean(selectedMap(1:(nTrigs-1)/2,:));
        trigRed = mean(selectedMap((end-(nTrigs-1)/2)+1:end,:));
        trigColorsAtt2 = [trigBlue; trigRed];
        trigColorsPA4 = [.52 .37 .75; .31 .74 .40; .27 .51 .84; 1.0 .57 .22];
        
        set(0,'defaultLineLineWidth',1)
        
        %% plot figs
        condColors.Att = trigColorsAtt2;
        condColors.PA = trigColorsPA4;
        xlims = [t(1) t(end)];
        condYlims.wAmps = [-.5 2.5];
        condDiffLims.wAmps = [-.8 .8];
        condYlims.hAmps = [-5 40];
        condDiffLims.hAmps = [-5 5];
        
        switch measure
            case 'wAmpsAtt'
                m = 'wAmps';
                p = 'Att';
            case 'hAmpsAtt'
                m = 'hAmps';
                p = 'Att';
            case 'wAmpsPA'
                m = 'wAmps';
                p = 'PA';
            case 'hAmpsPA'
                m = 'hAmps';
                p = 'PA';
            otherwise
                error('measure not recognized')
        end
        ylims = condYlims.(m);
        diffLims = condDiffLims.(m);
        colors = condColors.(p);
        
        %% indiv subjects
        figure
        set(gcf,'Position',[30 450 1200 250])
        for iSubject = 1:nSubjects
            subplot(nSubjects,1,iSubject)
            hold on
            for i=1:size(valsMean,1)
                plot(t, vals(i,:,iSubject), 'color', colors(i,:))
            end
            plot(xlims, [0 0], 'k')
            xlim(xlims)
            ylim(ylims)
            for iEv = 1:numel(evTimes)
                vline(evTimes(iEv),'color','k','LineStyle',':');
            end
            if iSubject==1
                xlabel('time (ms)')
                ylabel('amplitude')
            end
        end
        
        if strcmp(p,'Att')
            figure
            set(gcf,'Position',[30 450 1200 250])
            for iSubject = 1:nSubjects
                subplot(1,nSubjects,iSubject)
                hold on
                plot(t, valsDiff(1,:,iSubject), 'k', 'LineWidth', 2)
                plot(xlims, [0 0], 'k')
                xlim(xlims)
                ylim(diffLims)
                for iEv = 1:numel(evTimes)
                    vline(evTimes(iEv),'color','k','LineStyle',':');
                end
                if iSubject==1
                    xlabel('time (ms)')
                    ylabel('amplitude difference (T2-T1)')
                end
            end
        end
        
        %% group
        figure
        hold on
        for i=1:size(valsMean,1)
            shadedErrorBar(t, valsMean(i,:), valsSte(i,:), {'color', colors(i,:), 'LineWidth', 3}, 1)
        end
        plot(xlims, [0 0], 'k')
        for iEv = 1:numel(evTimes)
            vline(evTimes(iEv),'color','k','LineStyle',':');
        end
        xlabel('time (ms)')
        ylabel('amplitude')
        
        if strcmp(p,'Att')
            figure
            hold on
            for i=1:2
                shadedErrorBar(t, valsMean(i,:), valsDiffSte, {'color', colors(i,:), 'LineWidth', 3}, 1)
            end
            plot(xlims, [0 0], 'k')
            for iEv = 1:numel(evTimes)
                vline(evTimes(iEv),'color','k','LineStyle',':');
            end
            xlabel('time (ms)')
            ylabel('amplitude')
            
            figure
            hold on
            shadedErrorBar(t, valsDiffMean, valsDiffSte, 'k', 1)
            plot(xlims, [0 0], 'k')
            ylim(diffLims)
            for iEv = 1:numel(evTimes)
                vline(evTimes(iEv),'color','k','LineStyle',':');
            end
            xlabel('time (ms)')
            ylabel('amplitude difference (T2-T1)')
        end
        
    case 'groupData'
        %% calculate group mean and ste
        fieldNames = fieldnames(groupData);
        nFields = numel(fieldNames);
        for iF = 1:nFields
            fieldName = fieldNames{iF};
            vals = groupData.(fieldName);
            sdim = numel(size(vals)); % subject dimension
            groupMean.(fieldName) = mean(vals, sdim);
            groupSte.(fieldName) = std(vals, 0, sdim)./sqrt(nSubjects);
        end
        
        %% plot figs
        switch measure
            case {'w','h'}
                rd_plotTADetectDiscrimGroupTS(A, measure, ...
                    groupData, groupMean, groupSte)
            case {'tf','stf'}
                rd_plotTADetectDiscrimGroupTimeFreq(A, measure, ...
                    groupMean.amps, groupMean.ampsAtt, groupMean.ampsPA, groupMean.paDiff)
            otherwise
                error('measure not recognized')
        end
    otherwise
        error('dataType not recognized')
end

% rd_timeFreqAnalysis.m

% single-trial target data
%% setup
exptDir = '/Local/Users/denison/Data/TAPilot/MEG';
sessionDir = 'R0890_20140806';
fileBase = 'R0890_TAPilot_8.06.14';
analStr = 'eti';

dataDir = sprintf('%s/%s', exptDir, sessionDir);

switch analStr
    case ''
        filename = sprintf('%s/%s.sqd', dataDir, fileBase);
        figDir = sprintf('%s/figures/raw', dataDir);
    otherwise
        filename = sprintf('%s/%s_%s.sqd', dataDir, fileBase, analStr);
        figDir = sprintf('%s/figures/%s', dataDir, analStr);
end
if ~exist(figDir,'dir')
    mkdir(figDir)
end

channelsL = [26 60 14 92]; % R0890
channelsR = [1 50 7 8]; % R0890
% channelsL = [92 60 15 14]; % R0817
% channelsR = [51 1 50 39]; % R0817

saveFigs = 1;

%% prepare trl
cfg = [];
cfg.dataset             = filename;
cfg.trialdef.prestim    = 0.5;
cfg.trialdef.poststim   = 1;
cfg.trialdef.trig       = 165:166;

nEventsExpected = 272;

[trl,Events] = mytrialfun_all(cfg,2.5,nEventsExpected);

%% collect trigger info
triggers_target = [Events.trigger]';

type = {Events.channel};
type = (cellfun(@str2num, type))';

trig_ind = (1:numel(triggers_target))';
trigger_info = [trig_ind,triggers_target,trl,type];

%% trial types
leftTrials = find(trigger_info(:,6)==165);
rightTrials = find(trigger_info(:,6)==166);

%% ft_preprocessing
cfg = [];
cfg.dataset     = filename;
cfg.channel     = 'MEG';
cfg.continuous  = 'yes';
cfg.trl         = trl;

data = ft_preprocessing(cfg);

%% ft_freqanalysis
% cfg.trials         = trials;
cfg.output         = 'pow';
cfg.method         = 'mtmconvol';
cfg.taper          = 'hanning';
cfg.foi            = 1:60;
cfg.t_ftimwin      = 10 ./ cfg.foi;
cfg.toi            = -0.5:0.01:1;
freq               = ft_freqanalysis(cfg, data);

%% plot one channel
chan = 14;
figure
imagesc(squeeze(freq.powspctrm(chan,:,:)))
set(gca,'XTick',1:25:numel(freq.time))
set(gca,'XTickLabel',freq.time(1:25:numel(freq.time)))
set(gca,'YDir','normal')
title(sprintf('channel %d', chan))
colorbar

%% average across channels
powspctrmL = squeeze(mean(freq.powspctrm(channelsL,:,:),1));
powspctrmR = squeeze(mean(freq.powspctrm(channelsR,:,:),1));

%% plot channel average
figure
imagesc(powspctrmL)
set(gca,'XTick',1:25:numel(freq.time))
set(gca,'XTickLabel',freq.time(1:25:numel(freq.time)))
set(gca,'YDir','normal')
title(['left channels ' num2str(channelsL)])
colorbar
if saveFigs
    rd_saveAllFigs(gcf,{'targetTimeFreqPickedChannelsL'},'im',figDir)
end

figure
imagesc(powspctrmR)
set(gca,'XTick',1:25:numel(freq.time))
set(gca,'XTickLabel',freq.time(1:25:numel(freq.time)))
set(gca,'YDir','normal')
title(['right channels ' num2str(channelsR)])
colorbar
if saveFigs
    rd_saveAllFigs(gcf,{'targetTimeFreqPickedChannelsR'},'im',figDir)
end

figure
idx0 = find(freq.time==0);
plot([powspctrmL(:,idx0) powspctrmR(:,idx0)])
legend('left channels','right channels')


ssvefFreqs = [15 20 30 40];
for iF = 1:numel(ssvefFreqs)
    idxSSVEF(iF) = find(cfg.foi==ssvefFreqs(iF));
end

figure
plot(freq.time, powspctrmL(idxSSVEF,:))
legend(num2str(ssvefFreqs'))
xlabel('time (s)')
title('left channels')

figure
plot(freq.time, powspctrmR(idxSSVEF,:))
legend(num2str(ssvefFreqs'))
xlabel('time (s)')
title('right channels')

% donders_playground.m

%% load data
load('/Volumes/RACHO/Hypatia_20160416/TADetectDiscrim/MEG/R0983_20150813/prep/R0983_TADeDi_8.13.15_ebi_prepCleanData.mat')

data = cleanPrepData;
clear cleanPrepData

cond = data.trial_info(:,4);

load data/data_hdr
cfg = [];
layout = ft_prepare_layout(cfg, data_hdr);

%% demean
cfg = [];
cfg.channel = {'MEG'};
cfg.demean = 'yes';
cfg.baselinewindow = [-0.2 0];

datab = ft_preprocessing(cfg, data);

%% ERF baselined data
cfg = [];
cfg.trials = find(cond~=167);
% cfg.trials = find(cond==167);

avgb = ft_timelockanalysis(cfg, datab);

%% multiplot baselined
cfg = [];
cfg.showlabels = 'yes';
cfg.layout = layout;
figure
ft_multiplotER(cfg, avgb) % use cell array to plot multiple conditions on same plot

%% ICA
%% resample data
cfg = [];
cfg.resamplefs = 150;
cfg.detrend = 'no';
data_downsamp = ft_resampledata(cfg,data);

%% run ICA
cfg = [];
cfg.method = 'runica';
comp = ft_componentanalysis(cfg, data_downsamp);

%% topoplot ICA
cfg = [];
cfg.component = 1:20;
cfg.layout = layout;
cfg.comment = 'no';
figure
ft_topoplotIC(cfg, comp)

%% look at time series
cfg = [];
cfg.channel = 1:10; % components to be plotted
cfg.viewmode = 'component';
cfg.layout = layout;
figure
ft_databrowser(cfg, comp)

%% ERF of ICA
cfg = [];
cfg.trials = find(cond~=167);

avg_comp = ft_timelockanalysis(cfg, comp);

%% look at average time series
cfg = [];
cfg.channel = 11:20; % components to be plotted
cfg.blocksize = 3.61; 
cfg.viewmode = 'component';
cfg.layout = layout;
ft_databrowser(cfg, avg_comp);

%% using icabrowser
% cfg = [];
% cfg.channel = 1:10;
% cfg.layout = layout;
% ft_icabrowser(cfg, comp);

%% time-freq of ICA
cfg = [];
cfg.output = 'fourier';
cfg.method = 'fft'; %mtmfft
cfg.taper = 'hanning'; %dpss
% cfg.tapsmofrq = 1; %%%% what is rpttap?

freq_avg_comp = ft_freqanalysis(cfg, avg_comp);
freq_comp = ft_freqanalysis(cfg, comp);

%% look at average spectra
% freq_comp2 = freq_comp;
% freq_comp2.dimord = 'chan_time';
% freq_comp2.fourierspctrm = squeeze(freq_comp.fourierspctrm);
% freq_comp2.time = freq_comp2.freq;

cfg = [];
cfg.channel = 1:10; % components to be plotted
cfg.blocksize = 3.61; 
cfg.viewmode = 'component';
cfg.layout = layout;
ft_databrowser(cfg, freq_comp); %%%% getting an error

figure
components = 1:5;
for i=1:numel(components)
    subplot(numel(components),1,i)
    plot(freq_comp.freq, abs(squeeze(mean(freq_comp.fourierspctrm(:,components(i),:)))));
    ylim([0 0.6e-16])
    enddss

figure
plot(freq_comp.freq, abs(squeeze(mean(freq_avg_comp.fourierspctrm(:,1:20,:),1))));









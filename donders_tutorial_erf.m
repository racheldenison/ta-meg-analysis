% tutorial_preproc

%% define trials
cfg = [];
cfg.dataset             = 'R1018_TADeDi_r10-14_11.18.15.sqd';
cfg.trialdef.prestim    = 0.5; % sec
cfg.trialdef.poststim   = 3.1;
cfg.trialdef.trig       = [161:164,167];
threshold               = 2.5;

% cfg = ft_definetrial(cfg);
[trl,events] = mytrialfun_all(cfg,threshold,[]);
cfg.trl      = trl;

% remove bad trials identified previously
% cfg.trl([2, 3, 4, 30, 39, 40, 41, 45, 46, 47, 51, 53, 59, 77, 85],:) = []; 

%% preprocess
cfg.channel = {'MEG'};
cfg.demean  = 'yes';
cfg.baselinewindow = [-0.2 0];

data = ft_preprocessing(cfg);

%% ERF
% calculate
cfg = [];
avg = ft_timelockanalysis(cfg, data);

% multiplot
load data/data_hdr.mat
cfg = [];
cfg.showlabels = 'yes'; 
cfg.fontsize = 6; 
cfg.layout = ft_prepare_layout(cfg, data_hdr);
cfg.ylim = [-3e-13 3e-13];
ft_multiplotER(cfg, avg); %%% what about for multiple conditions?

% topoplot
cfg = [];
cfg.xlim = [0.3 0.5];
cfg.colorbar = 'yes';
ft_topoplotER(cfg, avg);

% topoplot, multiple time points
cfg = [];
cfg.xlim = -0.2:0.1:1.0;  % Define 12 time intervals
cfg.zlim = [-2e-13 2e-13];      % Set the 'color' limits.
clf;
ft_topoplotER(cfg, avg);

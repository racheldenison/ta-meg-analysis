% multiplot_tester.m

%% example with my data
load data/data_hdr.mat
cfg = [];
layout = ft_prepare_layout(cfg, data_hdr);

foi = A.stfFoi;
cmap = flipud(lbmap(64,'RedBlue'));

% example data 1
toi = A.stfToi;
vals = groupMean.amps(:,:,1); % freq x time
sz = size(vals);
b = ones(1,sz(1),sz(2));
b(1,:,:) = vals;
powspctrm = repmat(b,157,1,1); % chan x freq x time

% example data 2
toi = A.stftwinvals;
vals = groupMean.PAAUT(:,:,1,1); % freq x time
sz = size(vals);
b = ones(1,sz(1),sz(2));
b(1,:,:) = vals;
powspctrm = repmat(b,157,1,1); % chan x freq x time

TFdata.label = data_hdr.label(1:157);
TFdata.dimord = 'chan_freq_time';
TFdata.freq = foi;
TFdata.time = toi;
TFdata.powspctrm = powspctrm;

cfg = [];
cfg.layout = layout;
cfg.colormap = cmap;
cfg.zlim = [-.06 .06];	

figure 
ft_multiplotTFR(cfg, TFdata);


%% ft example for reference
load ~/Desktop/dataFIC.mat
cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'MEG';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = 2:2:30;                         % analysis 2 to 30 Hz in steps of 2 Hz 
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.toi          = -0.5:0.05:1.5;                  % time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
TFRhann = ft_freqanalysis(cfg, dataFIC);

cfg = [];
cfg.baseline     = [-0.5 -0.1]; 
cfg.baselinetype = 'absolute'; 
cfg.zlim         = [-3e-27 3e-27];	        
cfg.showlabels   = 'yes';	
cfg.layout       = 'CTF151.lay';
figure 
ft_multiplotTFR(cfg, TFRhann);


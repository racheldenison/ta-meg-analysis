% wavelet_tests.m

%% average across trials
cfg = [];
cfg.removemean = 'no';
[timelock] = ft_timelockanalysis(cfg, prep_data);

%% make data from prep data average
data = prep_data;
data.time = [];
data.trial = [];

data.time{1} = timelock.time;
data.trial{1} = timelock.avg;

%% cfg
cfg = [];
cfg.channel    = 'MEG';	                
cfg.method     = 'wavelet';                
cfg.width      = 12; 
cfg.output     = 'pow';	
cfg.foi        = 30;	                
cfg.toi        = -0.5:0.001:3.1;	

channel = 15;

%% wavelet of all trials
wav1 = ft_freqanalysis(cfg, prep_data);

%% wavelet of average across trials
wav2 = ft_freqanalysis(cfg, data);

%% specest of average across trials
[spectrum,freqoi,timeoi] = ft_specest_wavelet(data.trial{1}, data.time{1}, ...
    'freqoi', cfg.foi, 'timeoi', cfg.toi, 'width', cfg.width);
pow = abs(spectrum).^2;

%% hilbert
Fbp = cfg.foi + [-1.6 1.6]; 
filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample, Fbp);
hil1 = ft_preproc_hilbert(filt).^2;
hil2 = abs(hilbert(filt(channel,:))).^2;

%% plots
nr = 5;
figure
subplot(nr,1,1)
plot(wav1.time, squeeze(wav1.powspctrm(channel,1,:)))
title('freqanalysis - all trials')

subplot(nr,1,2)
plot(wav2.time, squeeze(wav2.powspctrm(channel,1,:)))
title('freqanalysis - average across trials')

subplot(nr,1,3)
plot(timeoi, squeeze(pow(channel,1,:)))
title('specest - average across trials')

subplot(nr,1,4)
plot(data.time{1}, hil1(channel,:))
title('preproc_hilbert - average across trials')

subplot(nr,1,5)
plot(data.time{1}, hil2)
title('matlab hilbert - average across trials')

%% testing different freq bands for hilbert
figure
hold all
for i = 1:5;
    Fbp = cfg.foi + [-1.6 1.6]*i;
    filt = ft_preproc_bandpassfilter(data.trial{1}, data.fsample, Fbp);
    hiltest = abs(hilbert(filt(channel,:))).^2;
    plot(hiltest)
end
legend(num2str((1:i)'))

%% phase angle
ph = angle(squeeze(spectrum(channel,:,:)));
phc = phase(squeeze(spectrum(channel,:,:))); % continuous (straight line if accurate)

figure
subplot(2,1,1)
plot(data.time{1}, data.trial{1}(channel,:))
subplot(2,1,2)
plot(timeoi, ph)

%% check with previously saved condition data
% first load A
t = A.t/1000;
foi = 30;
width = 12;
response = squeeze(A.trigMean)';
[spectrum,freqoi,timeoi] = ft_specest_wavelet(response, t, 'freqoi', foi, 'width', width);
spec = squeeze(spectrum);

ph = angle(spec);
phc = [];
for iCond = 1:size(spec,1)
    phc(iCond,:) = phase(spec(iCond,:)); 
end

figure
subplot(3,1,1)
plot(timeoi, ph(1:end-1,:)')
title('stim')
subplot(3,1,2)
plot(timeoi, ph(end,:)','k')
title('blank')
subplot(3,1,3)
hold on
plot(timeoi, ph(1:end-1,:)')
plot(timeoi, ph(end,:)','k')
title('stim & blank')

figure
subplot(1,3,1)
plot(timeoi, phc(1:end-1,:)')
xlim(xlims)
ylim(ylims)
title('stim')
subplot(1,3,2)
plot(timeoi, phc(end,:)','k')
xlim(xlims)
ylim(ylims)
title('blank')
subplot(1,3,3)
hold on
plot(timeoi, phc(1:end-1,:)')
plot(timeoi, phc(end,:)','k')
xlim(xlims)
ylim(ylims)
title('stim & blank')




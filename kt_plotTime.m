function img = kt_plotTime(sqdfile, analStr, rejectTrials, rejectChannels)
% KT_PLOTTIME  Plots: (1) concatenated trial time series, (2) trial time
% series averaged by channel, (3) fft 
% 
% KT_PLOTTIME(sqdfile, analStr, rejectTrials, rejectChannels)
% 
% sqdfile = 'datafile.sqd'
% analStr = 'string of analysis steps'
% rejectTrials = 1 reject, 0 do not reject 
% rejectChannels = 1 reject, 0 do not reject 

%% setup
% session = 'R1507_20190hel627';
% session2= 'R1507_TA2_6.27.19';

% figDir = '/Users/kantian/Dropbox/Carrasco_Lab/ta_meg_kt/figures';
% plain 
% load(sprintf('/Users/kantian/Dropbox/Data/TA2/MEG/%s/prep/%s_ebi_prepCleanData.mat',session,session2)); 
% load(sprintf('/Users/kantian/Dropbox/Data/TA2/MEG/%s/prep/channels_rejected.mat',session))
% load(sprintf('/Users/kantian/Dropbox/Data/TA2/MEG/%s/prep/trials_rejected.mat',session))

% noise 
% load(sprintf('/Users/kantian/Dropbox/Data/TANoise/MEG/%s/prep/%s_ebi_prepCleanData.mat',session,session2)); 
% load(sprintf('/Users/kantian/Dropbox/Data/TANoise/MEG/%s/prep/channels_rejected.mat',session))
% load(sprintf('/Users/kantian/Dropbox/Data/TANoise/MEG/%s/prep/trials_rejected.mat',session))

%% read data 

dat = ft_read_data(sqdfile);
hdr = ft_read_header(sqdfile);

cfg                     = [];
cfg.dataset             = sqdfile;
cfg.trialdef.prestim    = 0.2; %0; %0.5; % sec
cfg.trialdef.poststim   = 2.3; %2.3; %3.1;
cfg.trialdef.trig       = [168,167]; %[168,167]; %[161:164,167]; %168 = precue, 167=blank
threshold = 2.5;
[trl,Events]            = mytrialfun_all(cfg,threshold,[]);

prep_data          = ft_preprocessing(struct('dataset',sqdfile,...
    'channel','MEG','continuous','yes','trl',trl));

input = prep_data; 

%% concatenate data 
dummy = [];
for iTrial = 1:length(input.trial) % 580 trials
        dummy2 = input.trial{1,iTrial}; % go through trial by trial cells 
        dummy = [dummy dummy2]; % concat horizontally into trial x channel matrix 
end
dataConcatTimeSeries = dummy; % use this as unfiltered concatenated data matrix Nchans x Ntime 

% channels 
nChannels = size(dataConcatTimeSeries,1); 
% trials 
nTrials = size(input.trial,2);

% time vectors 
tConcat = 1:length(dataConcatTimeSeries); 
tTrial = 1:2501; 

% plot to check 
% figure
% set(gcf, 'Position',  [100, 100, 1200, 500])
% plot(tConcat,dataConcatTimeSeries)
% title('concatenated trial time series')
% ylabel('time (frames, Fs = 1000)') 

%% low level hpf on concatenated data 

% toFilter = dataConcatTimeSeries; % NChans x NTime 
% Fsample = 1000; 
% Fhp = 0.2; % high pass frequency 
% N = 8259; % filter order (default 6) 
% type = 'firws'; 
% direc = 'onepass-zerophase';
% [concatFilt] = ft_preproc_highpassfilter(toFilter, Fsample, Fhp, N, type, direc); 
% 
% save(sprintf('%s_concatFilt.mat',session),'concatFilt')

%% reshape 

% reshape concat time series into trial time series 
toReshape = dataConcatTimeSeries'; 
reshape1 = toReshape;
reshape2 = reshape(reshape1,2501,nTrials,nChannels); 
reshape3 = permute(reshape2,[1 3 2]); 
reshape4 = mean(reshape3,3); 
trialTime = reshape3; % filtered time x channels x trials 

% save(sprintf('%s_trialTime_%s.mat', session, analstr),'trialTime')

%%  reject trials

if rejectTrials
trialTime(:,:,trials_rejected) = NaN; 
% trialTime(:,:,37) = NaN; 
end

%% reject channels

if rejectChannels
channels_rejected2 = char(channels_rejected);
channels_rejected3 = channels_rejected2(:,end-2:end); % cell to char, extract last three
channels_rejected4 = [];
for i = 1:size(channels_rejected3,1)
    channels_rejected4(i) = str2double([channels_rejected3(i,1) channels_rejected3(i,2) channels_rejected3(i,3)]); % to string
end
channels_rejected5 = channels_rejected4'; % so readable as index
trialTime(:,channels_rejected5,:) = NaN; 

rejectChTrialTime = permute(trialTime,[2,3,1]); % reshape to Channel x Trial x Time 

% rejectChTrialTime_ef = rejectChTrialTime; 
% save(sprintf('%s_rejectChTrialTime_%s.mat', session, analstr),'rejectChTrialTime_ef')

dataConcatTimeSeries(channels_rejected5,:) = NaN; 
end

%% do fft
Fsample = 1000; 
nfft = 4096; % nfft = 2^nextpow2; 
nT = nTrials; % number of trials 

% 
toFFTfilt = trialTime; 
% 
Y = fft(toFFTfilt,nfft)/nT; % scale by number of samples
f = Fsample/2*linspace(0,1,nfft/2+1); % Fs/2 is the maximum frequency that can be measured
amps = 2*abs(Y(1:nfft/2+1,:,:)); % Multiply by 2 since only half the energy is in the positive half of the spectrum?
ampsMean = nanmean(amps,3);

% plot to check 
figure
hold on
set(gcf, 'Position',  [100, 100, 800, 2000])

subplot 311 
plot(tConcat,dataConcatTimeSeries(:,:))
title('concatenated time series')
xlabel('Time (ms)')
ylabel('Amplitude')

subplot 312
trialavg = nanmean(trialTime,3); 
plot(tTrial,trialavg(:,:)) % trialed filtered data 
xlim([0,2501]); 
ylim([-6e-13 6e-13])
xlabel('Time (ms)')
ylabel('Amplitude')
title('trial time series')

subplot 313
loglog(f, ampsMean(:,:))
xlim([f(1) f(end)])
ylim([10e-15 10e-10])
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
title('FFT on concatenated data') 

img = gcf; 

% export_fig(sprintf('R1452_TA2_11.19.18_run02_bietfp_rejectCh'), '-png', '-p0.1', '-transparent');

%% avg by channel 
% 
% % avg by ch
% oneChMean = [];
% for iChannel = 1:nChannels %:154
%     oneCh = squeeze(rejectChTrialTime(iChannel,:,:)); % squeeze removes singleton dimensions   
%     oneChMean2 = nanmean(oneCh); 
%     oneChMean = [oneChMean; oneChMean2];
% %     ylim([-6e-13 4e-13])
% 
% end
% 
% figure % check averages by channel 
% imagesc(oneChMean)
% colorbar
% colormap('pink')
% title('avg ERP by channel')
% xlabel('time')
% ylabel('channel')
% export_fig(sprintf('%s_meanchannels_%s',session, analstr), '-png', '-p0.1', '-transparent');
% 
% %% finding peaks 
% 
% % T1 mean peak per channel
% window1 = 1200:1500; % window to check for first max 
% [peak1, idx1] = max(oneChMean(:,window1),[],2); % returns max within defined window, and index relative to window
% idx1 = idx1 + min(window1) - 1; % convert window to trial time 
% 
% % T2 mean peak per channel
% window2 = 1501:1800; % window to check for second max 
% [peak2, idx2] = max(oneChMean(:,window2),[],2); % returns max within defined window, and index relative to window
% idx2 = idx2 + min(window2) - 1; % convert window to trial time 
% 
% figure
% % plot peaks
% for iChannel = 1 % just 1 channel or many [125 138 90 123 106 130 154 139 43 151]
%     plot(tTrial,oneChMean(iChannel,:))
%     % ylim([-6e-13 4e-13])
%     hold on
%     set(gcf, 'Position',  [100, 100, 1000, 500])
% 
%     title(sprintf('Channel %d',iChannel))
%     xlabel('Time (ms)]')
%     ylabel('Amplitude')
%     
%     vline(idx1(iChannel),'r') % vertical line at max1
%     label1 = sprintf('peak 1 %d at %dms', peak1(iChannel), idx1(iChannel));
%     text(idx1(iChannel),0,label1) % -1.2e-13
%     
%     vline(idx2(iChannel),'r') % vertical line at max2
%     label2 = sprintf('peak 2 %d at %dms', peak2(iChannel), idx2(iChannel));
%     text(idx2(iChannel),1e-13,label2) %  -1e-13
%     
%     vline(200,'k'); % baseline window check
%     hold off
%     
%     pause(2)
%     % export_fig(sprintf('%s_channel%d', session, iChannel), '-png', '-p0.1', '-transparent');
% end
% export_fig(sprintf('peak_%s',analstr), '-png', '-p0.1', '-transparent');
% 
% % plot check peak timing 
% figure
% hold on
% binedges = tTrial; 
% histogram(idx1, binedges, 'FaceColor', [0 0.4470 0.7410], 'FaceAlpha', 0.5, 'EdgeColor', 'none')
% histogram(idx2, binedges, 'FaceColor', [0.8500 0.3250 0.0980], 'FaceAlpha', 0.5, 'EdgeColor', 'none')
% title('check T1 T2 peak timing')
% legend('peak T1', 'peak T2')
% xlabel('time')
% ylabel('count (nTrials)')
% export_fig(sprintf('%s_checkpeaks_%s', session, analstr), '-png', '-p0.1', '-transparent');
% 
% % add save idx peak time. keep peak time same throughout. 
% 
% %% find variance per channel at peak 1 time 
% 
% datpeak1 = []; % channel * trial 
% datpeak1Unfilt = [];
% for iChannel = 1:nChannels 
%     dummy = rejectChTrialTime(iChannel,:,idx1(iChannel)); % erp at peak time for all trials by channel 
%     datpeak1  = [datpeak1 ; dummy];
% end
% 
% stdpeak1 = nanstd(datpeak1,[],2);
% 
% snr1 = peak1./stdpeak1; 
% [snr1sort snr1sort_idx] = sort(snr1,'descend');
% 
% % save
% snr1_ef = snr1; 
% snr1sort_ef = snr1sort; 
% snr1sort_idx_ef = snr1sort_idx; 
% 
% save(sprintf('%s_snr1_%s.mat', session, analstr),'snr1_ef')
% save(sprintf('%s_snr1sort_%s.mat', session, analstr),'snr1sort_ef')
% save(sprintf('%s_snr1sort_idx_%s.mat', session, analstr),'snr1sort_idx_ef')
% 
% % plot to check 
% colors = distinguishable_colors(12,[1 1 1; 0 0 0]); % generate perceptually distinct colors 
% x = 1:nChannels; 
% 
% % scatterplots of snr filtered v unfiltered 
% figure
% hold on
% set(gcf, 'Position',  [100, 100, 1000, 400])
% 
% % 
% load('R1507_20190627_snr1sort_e.mat')
% load('R1507_20190627_snr1sort_ef.mat')
% load('R1507_20190627_snr1sort_efl.mat')
% load('R1507_20190627_snr1sort_eflb.mat')
% load('R1507_20190627_snr1sort_eflbi.mat')
% %
% scatter(x,snr1sort_e',20,'filled','MarkerFaceColor',colors(1,:),'MarkerFaceAlpha',.7) 
% scatter(x,snr1sort_ef',20,'filled','MarkerFaceColor',colors(2,:),'MarkerFaceAlpha',.7)
% scatter(x,snr1sort_efl',20,'filled','MarkerFaceColor',colors(3,:),'MarkerFaceAlpha',.7)
% scatter(x,snr1sort_eflb',20,'filled','MarkerFaceColor',colors(4,:),'MarkerFaceAlpha',.7)
% scatter(x,snr1sort_eflbi',20,'filled','MarkerFaceColor',colors(5,:),'MarkerFaceAlpha',.7)
% scatter(x,snr1sort_flbtic',20,'filled','MarkerFaceColor',colors(6,:),'MarkerFaceAlpha',.7)
% % 
% legend('e','ef','efl','eflb','eflbi','flbtic')
% title(sprintf('%s \n SNR',session))
% ylabel('mean/std at time of mean peak1') 
% xlabel('channels')
% export_fig(sprintf('%s_scatterSNR_%s',session,analstr), '-png', '-p0.1', '-transparent');
% 
% %% boxplot preprocessing stepwise snrs 
% snrs = [snr1sort_ef snr1sort_efl snr1sort_eflb snr1sort_eflbi];
% 
% figure
% set(gcf, 'Position',  [100, 100, 400, 600])
% hold on 
% boxplot(snrs)
% title(sprintf('%s \n SNR',session))
% ylabel('peak1/std at mean peak1 time') 
% set(gca, 'XTickLabel', {'ef','efl','eflb','eflbi'})  
% xtickangle(45)
% 
% export_fig(sprintf('%s_boxplotSNR_flbi',session), '-png', '-p0.1', '-transparent');
% 
% %% single trial comparison plots before/after preproc step
% 
% % if importing matrices post 
% % snr1sort_idx = snr1sort_idx_e; 
% % snr1sort = snr1sort_e; 
% 
% colors = distinguishable_colors(6,[1 1 1; 0 0 0]); % generate perceptually distinct colors 
% 
% load('R1507_20190627_rejectChTrialTime_e.mat')
% load('R1507_20190627_rejectChTrialTime_ef.mat')
% load('R1507_20190627_rejectChTrialTime_efl.mat')
% % load('R1507_20190627_rejectChTrialTime_eflb.mat')
% load('R1507_20190627_rejectChTrialTime_eflbi.mat')
% load('R1507_20190627_rejectChTrialTime_eflbic.mat')
% load('R1507_20190627_rejectChTrialTime_flbtic.mat')
% 
% trialTime_none = permute(rejectChTrialTime,[3,1,2]); 
% trialTime_e = permute(rejectChTrialTime_e,[3,1,2]); 
% trialTime_ef = permute(rejectChTrialTime_ef,[3,1,2]);
% trialTime_efl = permute(rejectChTrialTime_efl,[3,1,2]);
% % trialTime_eflb = permute(rejectChTrialTime_eflb,[3,1,2]);
% trialTime_eflbi = permute(rejectChTrialTime_eflbi,[3,1,2]);
% trialTime_eflbic = permute(rejectChTrialTime_eflbic,[3,1,2]);
% trialTime_flbtic = permute(rejectChTrialTime_flbtic,[3,1,2]);
% 
% pickedChannels = 111;% snr1sort_idx(10:13); 
% pickedTrials = 100; % [10,30,60,100,200,223,500]; %[10,30,50,60,100,200,223,300,400,500];
% 
% figure
% plotpos = 1; 
% set(gcf,'Position',[100 100 2400 1200])
% for iChannel = 1:numel(pickedChannels) 
%     for iTrial = 1:numel(pickedTrials)       
%     subplot (numel(pickedChannels),numel(pickedTrials),plotpos)
%     hold on
%     yline(0,'k');
%     xlim([2,2501])
%     
%     plot(trialTime_e(:,pickedChannels(iChannel),pickedTrials(iTrial)),'Color',colors(1,:))
%     plot(trialTime_ef(:,pickedChannels(iChannel),pickedTrials(iTrial)),'Color',colors(2,:))
%     plot(trialTime_efl(:,pickedChannels(iChannel),pickedTrials(iTrial)),'Color',colors(3,:))
%     % plot(trialTime_eflb(:,pickedChannels(iChannel),pickedTrials(iTrial)))
%     plot(trialTime_eflbi(:,pickedChannels(iChannel),pickedTrials(iTrial)),'Color',colors(4,:))
%     plot(trialTime_eflbic(:,pickedChannels(iChannel),pickedTrials(iTrial)),'Color',colors(5,:))
%     plot(trialTime_flbtic(:,pickedChannels(iChannel),pickedTrials(iTrial)),'Color',colors(6,:))
%     title(sprintf('channel %d trial %d',pickedChannels(iChannel),pickedTrials(iTrial)))
%     plotpos = plotpos + 1; 
%     end
% end
% 
% legend('zero','e','ef','efl','eflbi','eflbic','flbtic') % ('e','ef','efl','eflb','eflbi')
% export_fig(sprintf('%s_singletrialchannels', session), '-png', '-p0.1', '-transparent');
% 
% % selective 
% figure
% plotpos = 1; 
% set(gcf,'Position',[100 100 2400 1200])
% for iChannel = 1:numel(pickedChannels) 
%     for iTrial = 1:numel(pickedTrials)       
%     subplot (numel(pickedChannels),numel(pickedTrials),plotpos)
%     hold on
%     yline(0,'k');
%     xlim([2,2501])
%     
%     % plot(trialTime_e(:,pickedChannels(iChannel),pickedTrials(iTrial)),'Color',colors(1,:))
%     % plot(trialTime_ef(:,pickedChannels(iChannel),pickedTrials(iTrial)),'Color',colors(2,:))
%     % plot(trialTime_efl(:,pickedChannels(iChannel),pickedTrials(iTrial)),'Color',colors(3,:))
%     % plot(trialTime_eflb(:,pickedChannels(iChannel),pickedTrials(iTrial)))
%     plot(trialTime_eflbi(:,pickedChannels(iChannel),pickedTrials(iTrial)),'Color',colors(4,:))
%     pause(0.5)
%     plot(trialTime_eflbic(:,pickedChannels(iChannel),pickedTrials(iTrial)),'Color',colors(5,:))
%     pause(0.5)
%     plot(trialTime_flbtic(:,pickedChannels(iChannel),pickedTrials(iTrial)),'Color',colors(6,:))
%     pause(0.5)
%     title(sprintf('channel %d trial %d',pickedChannels(iChannel),pickedTrials(iTrial)))
%     plotpos = plotpos + 1; 
%     end
% end
% 
% %%  whole sequence 1 channel all trials 
% figure
% plotpos = 1; 
% npreProcSteps = 6; 
% set(gcf,'Position',[100 100 2400 1200])
% for iChannel = 111%:numel(pickedChannels)    
%     subplot (npreProcSteps,1,plotpos)
%     hold on
%     yline(0,'k');
%     xlim([2,2501])
%     toPlot = squeeze(trialTime_eflbic(:,iChannel,:));
%     plot(toPlot)
%     title(sprintf('channel %d eflbic',iChannel))
%     plotpos = plotpos + 1; 
%     
%     subplot (npreProcSteps,1,plotpos)
%     hold on
%     yline(0,'k');
%     xlim([2,2501])
%     toPlot = squeeze(trialTime_eflbi(:,iChannel,:));
%     plot(toPlot)
%     title(sprintf('channel %d eflbi',iChannel))
%     plotpos = plotpos + 1; 
%     
%     subplot (npreProcSteps,1,plotpos)
%     hold on
%     yline(0,'k');
%     xlim([2,2501])
%     toPlot = squeeze(trialTime_flbtic(:,iChannel,:));
%     plot(toPlot)
%     title(sprintf('channel %d flbtic',iChannel))
%     plotpos = plotpos + 1; 
%     
%     subplot (npreProcSteps,1,plotpos)
%     hold on
%     yline(0,'k');
%     xlim([2,2501])
%     toPlot = squeeze(trialTime_efl(:,iChannel,:));
%     plot(toPlot)
%     title(sprintf('channel %d efl',iChannel))
%     plotpos = plotpos + 1; 
%     
%     subplot (npreProcSteps,1,plotpos)
%     hold on
%     yline(0,'k');
%     xlim([2,2501])
%     toPlot = squeeze(trialTime_ef(:,iChannel,:));
%     plot(toPlot)
%     title(sprintf('channel %d ef',iChannel))
%     plotpos = plotpos + 1; 
%     
%     subplot (npreProcSteps,1,plotpos)
%     hold on
%     yline(0,'k');
%     xlim([2,2501])
%     toPlot = squeeze(trialTime_e(:,iChannel,:));
%     plot(toPlot)
%     title(sprintf('channel %d e',iChannel))
%     plotpos = plotpos + 1; 
% end
% 
% export_fig(sprintf('%s_singletrialchannels_Channel%dTrialAll', session,iChannel), '-png', '-p0.1', '-transparent');
% 
% %%  whole sequence 1 channel 1 trial 
% figure
% plotpos = 1; 
% npreProcSteps = 6; 
% set(gcf,'Position',[100 100 2400 1200])
% for iChannel = 111%:numel(pickedChannels)
%     for iTrial = 100
%     subplot (npreProcSteps,1,plotpos)
%     hold on
%     yline(0,'k');
%     xlim([2,2501])
%     toPlot = squeeze(trialTime_eflbic(:,iChannel,iTrial));
%     plot(toPlot)
%     title(sprintf('channel %d trial %d eflbic',iChannel,iTrial))
%     plotpos = plotpos + 1; 
%     
%     subplot (npreProcSteps,1,plotpos)
%     hold on
%     yline(0,'k');
%     xlim([2,2501])
%     toPlot = squeeze(trialTime_eflbi(:,iChannel,iTrial));
%     plot(toPlot)
%     title(sprintf('channel %d trial %d eflbi',iChannel,iTrial))
%     plotpos = plotpos + 1; 
%     
%     subplot (npreProcSteps,1,plotpos)
%     hold on
%     yline(0,'k');
%     xlim([2,2501])
%     toPlot = squeeze(trialTime_flbtic(:,iChannel,iTrial));
%     plot(toPlot)
%     title(sprintf('channel %d trial %d flbtic',iChannel,iTrial))
%     plotpos = plotpos + 1; 
%     
%     subplot (npreProcSteps,1,plotpos)
%     hold on
%     yline(0,'k');
%     xlim([2,2501])
%     toPlot = squeeze(trialTime_efl(:,iChannel,iTrial));
%     plot(toPlot)
%     title(sprintf('channel %d trial %d efl',iChannel,iTrial))
%     plotpos = plotpos + 1; 
%     
%     subplot (npreProcSteps,1,plotpos)
%     hold on
%     yline(0,'k');
%     xlim([2,2501])
%     toPlot = squeeze(trialTime_ef(:,iChannel,iTrial));
%     plot(toPlot)
%     title(sprintf('channel %d trial %d ef',iChannel,iTrial))
%     plotpos = plotpos + 1; 
%     
%     subplot (npreProcSteps,1,plotpos)
%     hold on
%     yline(0,'k');
%     xlim([2,2501])
%     toPlot = squeeze(trialTime_e(:,iChannel,iTrial));
%     plot(toPlot)
%     title(sprintf('channel %d trial %d e',iChannel,iTrial))
%     plotpos = plotpos + 1; 
%     end
% end
% 
% export_fig(sprintf('%s_singletrialchannels_Channel%dTrial%d', session,iChannel,iTrial), '-png', '-p0.1', '-transparent');


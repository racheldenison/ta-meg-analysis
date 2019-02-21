% rd_TA2.m

%% i/o
exptDir = '/Local/Users/denison/Data/TA2/MEG';
sessionDir = 'R1187_20181119';

dataFile = dir(sprintf('%s/%s/mat/*condData.mat', exptDir, sessionDir));
dataFileName = sprintf('%s/%s/mat/%s', exptDir, sessionDir, dataFile.name);
figDir = sprintf('%s/%s/figures/ebi_ft', exptDir, sessionDir);

saveFigs = 0;

%% load data
load(dataFileName)

% load data header for plotting topologies
load data/data_hdr.mat

%% setup
t = D.t;
Fs = D.Fs;
eventTimes = D.eventTimes;
data = D.condData;

nT = numel(t);
sz = size(data);
nCue = sz(1);
nT1 = sz(2);
nT2 = sz(3);

cueNames = {'precue T1','precue T2','neutral'};

%% high-pass filtered time series
samplingInterval = 1;
tau = 100;
filtTau = samplingInterval/tau;

% Fhp = 0.1;
% % N = 4;
% type = 'fir';
% % type = 'firws';
% % dir = 'onepass-zerophase';
% % df = .01;

dataRaw = [];
dataFilt = [];
for iCue = 1:nCue
    for iT1 = 1:nT1
        for iT2 = 1:nT2
            fprintf('.')
            nTrials = size(data{iCue,iT1,iT2},3);
            for iTrial = 1:nTrials
                vals = data{iCue,iT1,iT2}(:,:,iTrial);
                idx = isnan(vals(1,:));
                vals(:,idx) = [];
                
                dataRaw{iCue,iT1,iT2}(:,~idx,iTrial) = vals;
                
                % time constant method
                valsfilt = filter([1-filtTau filtTau-1],[1 filtTau-1], vals);
                dataFilt{iCue,iT1,iT2}(:,~idx,iTrial) = valsfilt;

%                 % ft method
% %                 valsfilt = ft_preproc_highpassfilter(vals',Fs,Fhp,N,type,dir,[],df);
%                 valsfilt = ft_preproc_highpassfilter(vals',Fs,Fhp,[],type);
%                 dataFilt{iCue,iT1,iT2}(:,~idx,iTrial) = valsfilt';
            end
        end
    end
end

%% band-pass filtered time series
% Fbp = [1 40];
% 
% dataFilt = [];
% for iCue = 1:nCue
%     for iT1 = 1:nT1
%         for iT2 = 1:nT2
%             nTrials = size(data{iCue,iT1,iT2},3);
%             for iTrial = 1:nTrials
%                 vals = data{iCue,iT1,iT2}(:,:,iTrial);
%                 idx = isnan(vals(1,:));
%                 if all(idx)
%                     dataFilt{iCue,iT1,iT2}(:,:,iTrial) = nan(size(vals));
%                 else
%                     vals(:,idx) = [];
%                     % valsfilt = bandpass(vals,Fbp,Fs);
%                     valsfilt = ft_preproc_bandpassfilter(vals,Fs,Fbp);
%                     dataFilt{iCue,iT1,iT2}(:,~idx,iTrial) = valsfilt;
%                 end
%             end
%         end
%     end
% end

%% trial-by-trial baseline correction
twin = [800 1000]; %[-200 0];
tidx = find(t==twin(1)):find(t==twin(2));

dataInput = dataFilt; % dataRaw, dataFilt;

dataB = [];
for iCue = 1:nCue
    for iT1 = 1:nT1
        for iT2 = 1:nT2
            vals = dataInput{iCue,iT1,iT2};
            baselineVals = mean(vals(tidx,:,:),1);
            dataB{iCue,iT1,iT2} = vals - repmat(baselineVals,nT,1,1);
        end
    end
end

%% average across trials
dataInput = dataB; % dataRaw, dataB, dataFilt

dataMean = [];
for iCue = 1:nCue
    for iT1 = 1:nT1
        for iT2 = 1:nT2
            vals = dataInput{iCue,iT1,iT2};
            dataMean(:,:,iCue,iT1,iT2) = nanmean(vals,3);
        end
    end
end

trigMean = dataMean(:,:,:);

%% average by cue
dataMeanAtt = mean(mean(dataMean,5),4);

%% average across conditions 
% note, different from averaging across all single trials
trigMeanAll = mean(trigMean,3);

twin = [1000 1700];
tidx = find(t==twin(1)):find(t==twin(2));
maxval = max(max(abs(trigMeanAll(tidx,:))));

%% plot all channels
figure
subplot(2,1,1)
hold on
plot(t, trigMeanAll)
xlim([0 2300])
ylim([-maxval maxval]*1.1)
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end

subplot(2,1,2)
hold on
plot(t, abs(trigMeanAll))
xlim([0 2300])
ylim([0 maxval]*1.1)
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
xlabel('time (ms)')
ylabel('amplitude')

if saveFigs
    rd_saveAllFigs(gcf, {'tsAbsTs'}, 'plot', figDir);
end

%% find peak
% grandAbsMean = mean(abs(trigMeanAll),2);

%% find visually responsive channels
thresh = 100; % data = 170, dataFB = 100 
twin = [1000 1700];
tidx = find(t==twin(1)):find(t==twin(2));

channels = find(any(abs(trigMeanAll(tidx,:))>thresh,1));

% plot
fH(1) = figure;
hold on
plot(t, trigMeanAll(:,channels))
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
xlim([0 2300])
xlabel('time (ms)')
ylabel('amplitude')

% plot
spacer = 100;
fH(2) = figure;
hold on
plot(t, trigMeanAll(:,channels) + repmat(spacer*(1:numel(channels)),nT,1))
for iEv = 1:numel(eventTimes)
    vline(eventTimes(iEv),'k');
end
xlim([0 2300])
xlabel('time (ms)')

if saveFigs
    rd_saveAllFigs(fH, {'ts_visChannels','ts_vsChannels_stacked'}, 'plot', figDir);
end

%% find visually responsive peak
grandAbsMeanChannels = mean(abs(trigMeanAll(:,channels)),2);
grandAbsMaxChannels = max(abs(trigMeanAll(:,channels)),[],2);

twin = [1000 1200];
tidx = find(t==twin(1)):find(t==twin(2));
peakVal = max(grandAbsMaxChannels(tidx));
peakIdx = find(grandAbsMaxChannels==peakVal);
peakTime = t(peakIdx);
lims = [-peakVal peakVal]*1.1;

%% topo plot
vals = trigMeanAll(t==peakTime,:);
clims = lims;

[~, highChannel] = max(vals);
[~, lowChannel] = min(vals);

[~, channelOrder] = sort(vals);
lowChannels = channelOrder(1:3);
highChannels = channelOrder(end:-1:end-2);

figure
fH = ssm_plotOnMesh(vals, '', [], data_hdr, '2d',[],'numbers');
set(gca,'CLim',clims)
colorbar
title(sprintf('%d ms', peakTime))

if saveFigs
    rd_saveAllFigs(gcf, {sprintf('%dms',peakTime)}, 'topo', figDir);
end

%% plot selected channels
channels = [highChannel lowChannel];
nCh = numel(channels);

fH(1) = figure;
for iCh = 1:nCh
    subplot(nCh,1,iCh)
    plot(t, squeeze(dataMeanAtt(:,channels(iCh),:)))
    xlim([1000 1800])
    ylim(lims)
    if iCh==1
        legend(cueNames)
    end
    if iCh~=nCh
        set(gca,'XTickLabel','')
        set(gca,'YTickLabel','')
    end
    title(sprintf('ch %d',channels(iCh)))
end
xlabel('time (ms)')
ylabel('amplitude')

fH(2) = figure;
subplot(2,1,1)
plot(t, squeeze(mean(dataMeanAtt(:,highChannels,:),2)))
xlim([1000 1800])
ylim(lims)
title('high channels')
legend(cueNames)
subplot(2,1,2)
plot(t, squeeze(mean(dataMeanAtt(:,lowChannels,:),2)))
xlim([1000 1800])
ylim(lims)
title('low channels')
xlabel('time (ms)')
ylabel('amplitude')

if saveFigs
    rd_saveAllFigs(fH, {'tsAttHighLow','tsAttHighLowAve'}, 'plot', figDir);
end

%% topo movie
% clims = [-150 150];
% ts = 900:10:1800;
% 
% figure
% for iT = 1:numel(ts)
%     vals = trigMeanAll(t==ts(iT),:);
%     ssm_plotOnMesh(vals, '', [], data_hdr, '2d');
%     set(gca,'CLim',clims)
%     colorbar
%     title(sprintf('t = %d', ts(iT)))
%     pause(0.1)
% end

%% topo movie T1 and T2
clims = [-150 150];
ts = [1100:10:1400; 1400:10:1700];

figure('Position',[250 850 950 450])
for iT = 1:size(ts,2)
    for iEv = 1:size(ts,1)
        vals = trigMeanAll(t==ts(iEv,iT),:);
        subplot(1,size(ts,1),iEv)
        ssm_plotOnMesh(vals, '', [], data_hdr, '2d');
        set(gca,'CLim',clims)
        colorbar
        title(sprintf('t = %d', ts(iEv,iT)))
    end
    pause(0.2)
%     input('go')
end

%% topo different conditions
attT1_T1V = [1 7];
attT2_T1V = [2 8];
attT1_T1H = [4 10];
attT2_T1H = [5 11];
condSets = {attT1_T1V, attT2_T1V, attT1_T1H, attT2_T1H};
clims = [-peakVal peakVal];

figure
for i = 1:numel(condSets)
    subplot(2,2,i)
    vals = mean(trigMean(t==peakTime,:,condSets{i}),3);
    ssm_plotOnMesh(vals, '', [], data_hdr, '2d');
    set(gca,'CLim',clims)
    if i==1
        title('precue T1')
        ylabel('V')
    elseif i==2
        title('precue T2')
    elseif i==3
        ylabel('H')
    end
end

%% decoding
dataInput = dataRaw; % dataRaw, dataFilt, dataB
target = 'T1';
targetWindows = {[1000 1400],[1300 1700]};
sp = 5; % sampling period
kfold = 5;
svmops = sprintf('-s 0 -t 0 -c 1 -v %d -q', kfold);

switch target
    case 'T1'
        twin = targetWindows{1};
        d1 = dataInput(:,1,:); % T1 vertical
        d2 = dataInput(:,2,:); % T1 horizontal
    case 'T2'
        twin = targetWindows{2};
        d1 = dataInput(:,:,1); % T2 vertical
        d2 = dataInput(:,:,2); % T2 horizontal
end
d1 = d1(:); d2 = d2(:);

times = twin(1):sp:twin(2);
% times = 1000:sp:1400;
% time = peakTime;        

vals1 = []; vals2 = [];
for i = 1:numel(d1)
    vals1 = cat(3, vals1, d1{i}); 
    vals2 = cat(3, vals2, d2{i}); 
end

vals0 = cat(3, vals1, vals2);
labels0 = [ones(size(vals1,3),1); zeros(size(vals2,3),1)];

% stratify
nSamples = numel(labels0);
foldSize = nSamples/kfold/2; % 2 classes
stratIdx = [];
for iFold = 1:kfold
    idx1 = (1:foldSize) + (iFold-1)*foldSize;
    idx2 = idx1 + nSamples/2;
    stratIdx = [stratIdx idx1 idx2];
end

vals = vals0(:,:,stratIdx);
labels = labels0(stratIdx);
% vals = vals0;
% labels = labels0;


% classify
tic
acc = [];
for iTime = 1:numel(times)
    if mod(iTime,10)==0
        fprintf('%d, %s\n',iTime, datestr(now, 13))
    end
    
    time = times(iTime);
    % classification data
%     X = squeeze(vals(t==time,:,:))';
    X = squeeze(mean(vals(find(t==time):find(t==time+sp),:,:),1))'; % average across time window
    Y = labels;
    
    % remove nan
    idx = isnan(X(:,1));
    X(idx,:) = [];
    Y(idx) = [];
    
    % scale data
    Xs = zscore(X);
%     Xs = (X - repmat(min(X),size(X,1),1))./repmat(range(X),size(X,1),1); % [0 1]
    
    % fit and cross validate classifier
    acc(iTime) = svmtrain(Y, Xs, svmops); 
end
toc

% plot
xlims = times([1 end]);
ylims = [40 80];

figure
hold on
plot(times,acc)
plot(xlims,[50 50],'k')
xlim(xlims)
ylim(ylims)
xlabel('time (ms)')
ylabel('classification accuracy (%)')

% get the parameters of a sample model
model = svmtrain(Y,Xs,'-s 0 -t 0 -c 1');
w = model.SVs' * model.sv_coef;
b = -model.rho;
if (model.Label(1) == -1)
    w = -w; b = -b;
end

figure
fH = ssm_plotOnMesh(w', '', [], data_hdr, '2d',[],'numbers');
% set(gca,'CLim',clims)
colorbar






% rd_waveletAnalysis.m

%% load data
load /Local/Users/denison/Data/TAPilot/MEG/R0817_20140820/mat/R0817_TAPilot_8.20.14_erf_workspace.mat

channelsL = [92 60 15 14]; % R0817
channelsR = [51 1 50 39]; % R0817

%% wavelet
channels = channelsL;
for iChan = 1:numel(channels)
    data = squeeze(trigData(:,channels(iChan),:))'; % trials by samples
    [spectrum,freqoi,timeoi] = ft_specest_wavelet(data, t);
    specAmp = abs(squeeze(spectrum));
    
    freqIdx = find(abs(freqoi-ssvefFreq) == min((abs(freqoi-ssvefFreq))));
    
    wAmp = specAmp(freqIdx,:);
    wAmpNorm = wAmp./mean(wAmp(500:1000));
    wAmps(iTrial,:) = wAmpNorm;
end
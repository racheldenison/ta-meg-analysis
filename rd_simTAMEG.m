% rd_simTAMEG.m

%% setup
sinewave = @(A,f,t,ph) A*sin(2*pi*f*t + ph);

Fs = 1000;
dur = 5;

ssvefFreq = 20;
ssvefAmp = 1;
Fbp = ssvefFreq + [-1.6 1.6];

erfT = 2.3;
erfWidth = .2;
erfAmp = 1;

EndoT = 2;
EndoWidth = 0.1;
EndoAmp = 1;

noiseStd = 1;

t = 0:1/Fs:dur;

nTrials = 50;

plotFigs = 0;

%% run simulated trials
for iTrial = 1:nTrials
    %% make time series
    ssvef = sinewave(ssvefAmp,ssvefFreq,t,0);
    erf = makeGaussian(t,erfT,erfWidth,erfAmp) - makeGaussian(t,erfT+0.5*erfWidth,erfWidth,erfAmp);
    attnGain = makeGaussian(t,EndoT,EndoWidth,EndoAmp) + ones(size(t));
    noise = noiseStd.*randn(size(t));
    
    response = (ssvef + erf).*attnGain + noise;
    responses(iTrial,:) = response;
    
    if plotFigs
        figure
        subplot(2,1,1)
        hold on
        plot(t, ssvef)
        plot(t, erf, 'g')
        plot(t, attnGain, 'r')
        subplot(2,1,2)
        plot(t, response, 'k')
    end
    
    %% filter
%     [filt] = ft_preproc_bandpassfilter(dat,Fs,Fbp,N,type,dir,instabilityfix,df,wintype,dev,plotfiltresp)
    filt = ft_preproc_bandpassfilter(response,Fs,Fbp);
    
    if plotFigs
        figure
        subplot(2,1,1)
        hold on
        plot(t,ssvef.*attnGain)
        plot(t, filt, 'r')
        subplot(2,1,2)
        hold on
        plot(t,ssvef.*attnGain)
        plot(t, filt, 'r')
        xlim([1.5 2.5])
    end
    
    %% Hilbert transform
    h = hilbert(filt);
    hAmp = abs(h);
%     hAmp = ft_preproc_hilbert(filt, 'abs');
    hAmps(iTrial,:) = hAmp;
    
    if plotFigs
        figure
        plot(t, hAmp)
    end
end

%% average then Hilbert
meanResponse = mean(responses); % mean of raw time series
meanResponseF = ft_preproc_bandpassfilter(meanResponse,Fs,Fbp);
meanResponseFH = abs(hilbert(meanResponseF));

%% Hilbert then average
meanH = mean(hAmps);

%% average then Wavelet
% [spectrum,freqoi,timeoi] = ft_specest_wavelet(dat, time, varargin)
[spectrum,freqoi,timeoi] = ft_specest_wavelet(meanResponse, t);
specAmp = abs(squeeze(spectrum));

freqIdx = find(abs(freqoi-ssvefFreq) == min((abs(freqoi-ssvefFreq))));

wavAmp = specAmp(freqIdx,:);
wavAmpNorm = wavAmp./mean(wavAmp(500:1000));

if plotFigs
    figure
    subplot(2,1,1)
    imagesc(flipud(specAmp))
    subplot(2,1,2)
    plot(t, wavAmpNorm)
end

%% plots
% figure
% hold on
% plot(t, ssvef.*attnGain, 'r')
% plot(t, meanResponse)
% legend('orig signal', 'mean response')

figure
hold on
plot(t, attnGain, 'g')
plot(t, meanResponseFH)
plot(t, meanH, 'r')
plot(t, wavAmpNorm, 'k')
legend('true attn gain','average then Hilbert','Hilbert then average','average then wavelet')

    
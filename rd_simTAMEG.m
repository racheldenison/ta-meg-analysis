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

noiseStd = 0;

t = 0:1/Fs:dur;

nTrials = 2;

plotFigs = 1;

%% run simulated trials
for iTrial = 1:nTrials
    %% make time series
    ssvef = sinewave(ssvefAmp,ssvefFreq,t,0);
    erf = makeGaussian(t,erfT,erfWidth,erfAmp) - makeGaussian(t,erfT+0.5*erfWidth,erfWidth,erfAmp);
    attnGain = makeGaussian(t,EndoT,EndoWidth,EndoAmp) + ones(size(t));
    noise = noiseStd.*randn(size(t));
    
%     % contrast modulation from target
%     targetOn = find(t==3):find(t==3.05);
%     ssvef(targetOn) = ssvef(targetOn)*0.5;
    
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
    
    %% wavelet
    foi = ssvefFreq;
    width = 12;
%     [spectrum,freqoi,timeoi] = ft_specest_wavelet(response, t, 'freqoi', foi, 'width', width);
    [spectrum,freqoi,timeoi] = ft_specest_wavelet(response, t);
    specAmp = abs(squeeze(spectrum));
    
    freqIdx = find(abs(freqoi-ssvefFreq) == min((abs(freqoi-ssvefFreq))));
    
    wAmp = specAmp(freqIdx,:);
    wAmpNorm = wAmp./mean(wAmp(500:1000));
    wAmps(iTrial,:) = wAmpNorm;
end

%% Hilbert then average
meanH = mean(hAmps);

%% average then Hilbert
meanResponse = mean(responses); % mean of raw time series
meanResponseF = ft_preproc_bandpassfilter(meanResponse,Fs,Fbp);
meanResponseFH = abs(hilbert(meanResponseF));

%% wavelet then average
meanW = mean(wAmps);

%% average then wavelet
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
plot(t, meanW, 'c')
legend('true attn gain','average then Hilbert','Hilbert then average','average then wavelet','wavelet then average')

    
%% investigation of wavelet
foi = ssvefFreq;
width = 12;
[spectrum,freqoi,timeoi] = ft_specest_wavelet(response, t, 'freqoi', foi, 'width', width);
spec = squeeze(spectrum);

% phase example
x = cos(pi/4*(0:100));
y = hilbert(x);
sigphase = angle(y);
plot(sigphase)

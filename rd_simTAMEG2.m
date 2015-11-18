% rd_simTAMEG2.m

%% setup
sinewave = @(A,f,t,ph) A*sin(2*pi*f*t + ph);

Fs = 1000;
dur = 5;

ssvefFreq = 30;
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

plotFigs = 1;

%% make time series
ssvef = sinewave(ssvefAmp,ssvefFreq,t,0);
erf = makeGaussian(t,erfT,erfWidth,erfAmp) - makeGaussian(t,erfT+0.5*erfWidth,erfWidth,erfAmp);
attnGain = makeGaussian(t,EndoT,EndoWidth,EndoAmp) + ones(size(t));
noise = noiseStd.*randn(size(t));

%     % contrast modulation from target
%     targetOn = find(t==3):find(t==3.05);
%     ssvef(targetOn) = ssvef(targetOn)*0.5;

response = (ssvef + erf).*attnGain + noise;

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

%% Hilbert
filt = ft_preproc_bandpassfilter(response,Fs,Fbp);
hAmp = abs(hilbert(filt));

%% wavelet
foi = ssvefFreq;
width = 12;
[spectrum,freqoi,timeoi] = ft_specest_wavelet(response, t, 'freqoi', foi, 'width', width);
wAmp = abs(squeeze(spectrum));
wAmpNorm = wAmp./mean(wAmp(500:1000));
wAmp(500)

%% plots
figure
hold on
plot(t, attnGain, 'g')
plot(t, hAmp, 'r')
plot(t, wAmp, 'k')
plot(t, wAmpNorm, 'b')
legend('true attn gain','Hilbert','wavelet','wavelet normalized')

%% investigation of wavelet
% phase example
x = cos(pi/4*(0:100));
y = hilbert(x);
sigphase = angle(y);
% plot(sigphase)



% rd_simSampledSquarewaveConvBiphasic2.m

stimFreq = 20;

ks = 500:-50:50;
ns = 1:10;

%% make stim
[stim, t, Fs] = rd_simSampledSquarewaveStim(stimFreq);
stim = (stim+1)./2;
nSamples = numel(t);

figure
plot(t,stim,'.-')
xlim([0 1])

stimOn = stim;
stimOff = 1-stim;

figure
hold on
plot(t,stimOn)
plot(t,stimOff,'r')
xlim([0 1])

%% make biphasic impulse response function
for i = 1:numel(ks)
    for j = 1:numel(ns)
        k = ks(i);
        n = ns(j);
        g = makeBiphasic(0:.001:0.2,k,n);
%         g(g<0)=g(g<0)*0.5; % unbalanced
        g = g./max(g);
        gs(:,i,j) = g;
    end
end

figure
for i=1:numel(ks)
    subplot(numel(ns),1,i)
    plot(squeeze(gs(:,i,:)))
    xlim([0 200])
    set(gca,'YTickLabel','')
    if i~=numel(ks)
        set(gca,'XTickLabel','')
    else
        xlabel('time (ms)')
    end
end

%% convolve with impulse response function
for i = 1:numel(ks)
    for j = 1:numel(ns)
        g = gs(:,i,j);
        onResponses(:,i,j) = conv(stimOn, g);
        offResponses(:,i,j) = conv(stimOff, g);
    end
end
onResponses = onResponses(1:nSamples,:,:);
offResponses = offResponses(1:nSamples,:,:);

% rectify
onResponses(onResponses<0) = 0;
offResponses(offResponses<0) = 0;

% sum
responses = onResponses + offResponses;

%% example time series
figure
plot(t, responses(:,8,8));
xlim([0 1])
title('example')
xlabel('time (s)')

%% FFT
nfft = nSamples;
Y = fft(responses,nfft)/nSamples; % Scale by number of samples
f = Fs/2*linspace(0,1,nfft/2+1); % Fs/2 is the maximum frequency that can be measured
amps = 2*abs(Y(1:nfft/2+1,:,:)); % Multiply by 2 since only half the energy is in the positive half of the spectrum?

%% example spectrum
figure
plot(f, amps(:,1,5));
xlim([0 150])
title('example')

%% power at 1F and 2F
[xx, freqIdx(1)] = min(abs(f-stimFreq));
[xx, freqIdx(2)] = min(abs(f-stimFreq*2));

stimFreqAmps(:,:,1) = squeeze(amps(freqIdx(1),:,:)); % ks x ns
stimFreqAmps(:,:,2) = squeeze(amps(freqIdx(2),:,:));

figure
for iF = 1:numel(freqIdx)
    subplot(1,numel(freqIdx),iF)
    imagesc(stimFreqAmps(:,:,iF))
    clims = [0 8];
    set(gca,'clim',clims);
    colorbar
    xlabel('n parameter')
    ylabel('k parameter')
    title(sprintf('%d Hz', stimFreq*iF))
end

%% visualize gammas with high stim freq power
figure
for i=1:numel(ks)
    subplot(numel(ks),1,i)
    hold on
    for j=1:numel(ns)
        plot(gs(:,i,j),'color',1-repmat(stimFreqAmps(i,j,1)./3,1,3))
    end
    xlim([0 150])
    set(gca,'YTickLabel','')
    if i~=numel(ks)
        set(gca,'XTickLabel','')
    else
        xlabel('time (ms)')
    end
end


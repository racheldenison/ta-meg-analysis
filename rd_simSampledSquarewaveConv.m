% rd_simSampledSquarewaveConv.m

stimFreq = 20;

scales = 1:10;
shapes = 1:10;

%% make stim
[stim, t, Fs] = rd_simSampledSquarewaveStim(stimFreq);
stim = (stim+1)./2;
nSamples = numel(t);

figure
plot(t,stim,'.-')
xlim([0 1])

stimOn = double(diff(stim)==1);
stimOn(end+1) = 0;
stimOff = double(diff(1-stim)==1);
stimOff(end+1) = 0;

figure
hold on
plot(t,stimOn)
plot(t,stimOff,'r')
xlim([0 1])

%% make gammas
for i = 1:numel(scales)
    for j = 1:numel(shapes)
        scale = scales(i);
        shape = shapes(j);
        g = makeGamma(0:200, [], shape, scale, 1);
        gs(:,i,j) = g;
    end
end

figure
for j=1:numel(shapes)
    subplot(numel(shapes),1,j)
    plot(gs(:,:,j))
end

%% convolve with gamma
for i = 1:numel(scales)
    for j = 1:numel(shapes)
        g = gs(:,i,j);
        onResponses(:,i,j) = conv(stimOn, g);
        offResponses(:,i,j) = conv(stimOff, g);
    end
end
onResponses = onResponses(1:nSamples,:,:);
offResponses = offResponses(1:nSamples,:,:);
responses = onResponses + offResponses;

%% example time series
figure
plot(t, responses(:,3,5));
xlim([0 1])
title('example')
xlabel('time (s)')

%% FFT
nfft = 2^nextpow2(nSamples); % Next power of 2 from length of y
Y = fft(responses,nfft)/nSamples; % Scale by number of samples
f = Fs/2*linspace(0,1,nfft/2+1); % Fs/2 is the maximum frequency that can be measured
amps = 2*abs(Y(1:nfft/2+1,:,:)); % Multiply by 2 since only half the energy is in the positive half of the spectrum?

%% example spectrum
figure
plot(f, amps(:,1,1));
xlim([0 150])
title('example')

%% power at 1F and 2F
[xx, freqIdx(1)] = min(abs(f-stimFreq));
[xx, freqIdx(2)] = min(abs(f-stimFreq*2));

stimFreqAmps(:,:,1) = squeeze(amps(freqIdx(1),:,:)); % scales x shapes
stimFreqAmps(:,:,2) = squeeze(amps(freqIdx(2),:,:));

figure
for iF = 1:numel(freqIdx)
    subplot(1,numel(freqIdx),iF)
    imagesc(stimFreqAmps(:,:,iF))
%     clims = get(gca,'clim');
%     clims(1) = 0;
    clims = [0 0.5];
    set(gca,'clim',clims);
    colorbar
    xlabel('shape parameter')
    ylabel('scale parameter')
    title(sprintf('%d Hz', stimFreq*iF))
end

%% visualize gammas with high stim freq power
figure
for j=1:numel(shapes)
    subplot(numel(shapes),1,j)
    hold on
    for i=1:numel(scales)
        plot(gs(:,i,j),'color',1-repmat(stimFreqAmps(i,j,1)*2,1,3))
    end
    xlim([0 150])
    set(gca,'YTickLabel','')
    if j~=numel(shapes)
        set(gca,'XTickLabel','')
    else
        xlabel('time (ms)')
    end
end

% cutoff = max(max(stimFreqAmps(:,:,2)))/2; % set cutoff based on 2F
% goodGs = gs;
% excl = stimFreqAmps(:,:,1)<cutoff;
% goodGs(:,excl) = NaN;
% 
% figure
% for j=1:numel(shapes)
%     subplot(numel(shapes),1,j)
%     plot(goodGs(:,:,j))
%     xlim([0 150])
% end


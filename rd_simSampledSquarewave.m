% rd_simSampledSquarewave.m

%% setup
stimFreq = 20; % 15 or 20

dur = 2.5;
Fs = 120; % should be 120 if using the constructed neural response
t = 0:1/Fs:dur-1/Fs;
nSamples = numel(t);

%% constructed neural response - use 120 Hz sampling rate
% assumes that neurons respond to onsets and offsets
switch stimFreq
    case 15
        % one cycle of the stimulus at 60 Hz is [1 1 0 0]
        % so the neural response is [1 0 1 0]
        % sampled at 120 Hz, this would be [1 0 0 0 1 0 0 0]
        response = repmat([1 0 0 0],1,nSamples/4);
    case 20
        % one cycle of the stimulus at 60 Hz is [1 0 0]
        % so the neural response is [1 1 0]
        % sampled at 120 Hz, this would be [1 0 1 0 0 0]
        response = repmat([1 0 1 0 0 0],1,nSamples/6);
end

% figure
figure
plot(t, response, '.-')
xlabel('time (s)')
ylabel('response amp')
title(sprintf('%d Hz stimulus', stimFreq))
xlim([0 1])

%% FFT
nfft = 2^nextpow2(nSamples); % Next power of 2 from length of y
Y = fft(response,nfft)/nSamples; % Scale by number of samples
f = Fs/2*linspace(0,1,nfft/2+1); % Fs/2 is the maximum frequency that can be measured
amps = 2*abs(Y(1:nfft/2+1)); % Multiply by 2 since only half the energy is in the positive half of the spectrum?

% figure
figure
plot(f, amps)
xlabel('frequency')
ylabel('amplitude')
title(sprintf('%d Hz stimulus', stimFreq))

%% time-frequency on response
% taper          = 'hanning';
% foi            = 1:50;
% t_ftimwin      = 10 ./ foi;
% toi            = t;
% tfAmps = [];
% [spectrum,ntaper,freqoi,timeoi] = ft_specest_mtmconvol(response, t, ...
%     'timeoi', toi, 'freqoi', foi, 'timwin', t_ftimwin, ...
%     'taper', taper, 'dimord', 'chan_time_freqtap');
% specAmp = abs(squeeze(spectrum))';
% 
% % figure
% ytick = 10:10:numel(foi);
% xtick = 1:30:numel(toi);
% figure
% imagesc(specAmp)
% rd_timeFreqPlotLabels(toi,foi,xtick,ytick);

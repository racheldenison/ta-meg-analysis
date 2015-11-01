function [stim, t, Fs] = rd_simSampledSquarewaveStim(stimFreq)

%% setup
if nargin==0
    stimFreq = 20; % 15 or 20
end

dur = 2.4;
Fs = 1000;
t = 0:1/Fs:dur-1/Fs;
nSamples = numel(t);

plotFigs = 0;

%% stimulus
switch stimFreq
    case 20
        duty = 33;
    case 100
        duty = 10;
    otherwise
        duty = 50;
end
stim = square(2*pi*stimFreq*t, duty);

% figure
if plotFigs
    figure
    plot(t, stim, '.-')
    xlabel('time (s)')
    ylabel('stim amp')
    title(sprintf('%d Hz stimulus', stimFreq))
    xlim([0 1])
end

%% FFT
nfft = 2^nextpow2(nSamples); % Next power of 2 from length of y
Y = fft(stim,nfft)/nSamples; % Scale by number of samples
f = Fs/2*linspace(0,1,nfft/2+1); % Fs/2 is the maximum frequency that can be measured
amps = 2*abs(Y(1:nfft/2+1)); % Multiply by 2 since only half the energy is in the positive half of the spectrum?

% figure
if plotFigs
    figure
    plot(f, amps)
    xlabel('frequency')
    ylabel('amplitude')
    title(sprintf('%d Hz stimulus', stimFreq))
end

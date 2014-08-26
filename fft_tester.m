% fft_tester.m

plotFigs = 0;

Fs = 1000;
t = 0:1/Fs:50; % hmm, peak amp seems to depend on length of t
x = sin(2*pi*t); % sine wave frequency = 1

n = numel(t);
nfft = 2^nextpow2(n); % next power of 2 from length of x
Y = fft(x,nfft)/n; % scale by number of samples
f = Fs/2*linspace(0,1,nfft/2+1); % Fs/2 is the maximum frequency that can be measured
amps = 2*abs(Y(1:nfft/2+1)); % multiply by 2 since only half the energy is in the positive half of the spectrum?

if plotFigs
    figure
    subplot(2,1,1)
    plot(t,x)
    subplot(2,1,2)
    plot(f,amps)
end

peakAmp = max(amps)
peakFreq = f(amps==peakAmp)
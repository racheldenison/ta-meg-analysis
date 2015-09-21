% rd_simSampledSquarewavePhaseShift.m

%% setup
stimFreq = 20; % 15 or 20

dur = 2.4;
Fs = 1000; 
t = 0:1/Fs:dur-1/Fs;
nSamples = numel(t);

sinewave = @(A,f,t,ph) A*sin(2*pi*f*t + ph);

%% determine phase difference between on and off populations
switch stimFreq
    case 15
        duty = 0.5;
    case 20
        duty = 0.33;
end
phaseShift = 2*pi*duty;

%% make sine wave responses
onResponse = sinewave(1,stimFreq,t,0);
offResponse = sinewave(1,stimFreq,t,phaseShift);

% half-wave rectify
onResponse(onResponse < 0) = 0;
offResponse(offResponse < 0) = 0;

response = onResponse + offResponse;

%% FFT
nfft = 2^nextpow2(nSamples); % Next power of 2 from length of y
Y = fft(response,nfft)/nSamples; % Scale by number of samples
f = Fs/2*linspace(0,1,nfft/2+1); % Fs/2 is the maximum frequency that can be measured
amps = 2*abs(Y(1:nfft/2+1)); % Multiply by 2 since only half the energy is in the positive half of the spectrum?

%% figure
figure
subplot(3,1,1)
hold on
plot(t, onResponse)
plot(t, offResponse, 'r')
xlim([0 1])
ylim([0 2])
legend('on population','off population')
subplot(3,1,2)
plot(t, response, 'k')
xlim([0 1])
ylim([0 2])
legend('sum')
subplot(3,1,3)
plot(f, amps)
xlim([0 100])
xlabel('frequency')
ylabel('amplitude')



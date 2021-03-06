function [waveletDuration, spectralBandwidth] = rd_calculateWaveletResolution(F, width)

% [waveletDuration, spectralBandwidth] = rd_calculateWaveletResolution(F, width)
%
% F is the frequency in Hz
% width is the width of the wavelet in number of cycles
%
% one but not both inputs can be a vector

spectralBandwidth = F./width.*2;
waveletDuration = width./F./pi;
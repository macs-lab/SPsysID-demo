function [amplitude_ratio,phase_lag] = twoSine_amplitudeRatio_PhaseLag(dataIn,dataOut)
% function twoSine_amplitudeRatio_PhaseLag(dataIn,dataOut,Ts)
% computes the amplitude ratio and phase differnce between two sinusoidal
% signals at the same frequency
% inputs: 
%   dataIn: input sinusoidal signal
%   dataOut: output sinusoidal signal
%   Ts: sampling time
% outputs: 
%   amplitude_ratio: 
%           output magnitude 
%           ----------------
%           input magnitude
%   phase_lag: output phase - input phase
% 
% Xu Chen
% chx@uw.edu
% 2012-05-16
% 2015-02-08
% 2022-04-29

if nargin < 2
    error('insufficient inputs')
end

x = dataIn;
y = dataOut;

% remove bias
x = x - mean(x);
y = y - mean(y);

if length(x) ~= length(y)
    error('Lengths of inputs and outputs are not the same.')
end
npts = length(x);

% take the FFT
X=fft(x);
Y=fft(y);

% Calculate the numberof unique points
NumUniquePts = ceil((npts+1)/2);

% Determine the max value and max point.
% This is where the sinusoidal is located if the signal to noise ratio is large. 
[mag_x idx_x] = max(abs(X));
[mag_y idx_y] = max(abs(Y));

% determine the phase difference
% at the maximum point.
px = phase(X(idx_x));
py = phase(Y(idx_y));
phase_lag = py - px;
% phas_lag = mod(phase_lag+pi,2*pi)-pi;

% determine the amplitude scaling
amplitude_ratio = mag_y/mag_x;
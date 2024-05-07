function [amplitude_ratio,phase_lag] = twoSine_amplitudeRatio_PhaseLag_direct(dataIn,dataOut,freq_rad)
% function twoSine_amplitudeRatio_PhaseLag_direct(dataIn,dataOut,freq_rad)
% computes the amplitude ratio and phase differnce between two sinusoidal
% signals at the same frequency
% inputs: 
%   dataIn: input sinusoidal signal
%   dataOut: output sinusoidal signal
%   freq_rad: frequency of the sinusoidals in rad/sec
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

if nargin < 3
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

indx = 1:npts;
% basis vectors
sinb = sin(freq_rad*indx);
cosb = cos(freq_rad*indx);

sinb = reshape(sinb,size(x));
cosb = reshape(cosb,size(x));

reX = sum(x.*sinb);
imX = sum(x.*cosb);
AX = sqrt(reX^2+imX^2);
phaX = atan(imX/reX);

reY = sum(y.*sinb);
imY = sum(y.*cosb);
AY = sqrt(reY^2+imY^2);
phaY = atan(imY/reY);

amplitude_ratio = AY/AX;
phase_lag = phaY-phaX;

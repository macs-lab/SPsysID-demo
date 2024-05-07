function [gain,phase] = sysid_T2P(respT,freq_Hz,C,Ts)
% function [gain,phase] = sysid_T2P(respT,freq_Hz,C,Ts)
% transforms the frequency response data of a complementary sensitivity
% function T=PC/(1+PC) into the frequency response of the plant P
% required inputs:
%   respT: frequency response of the complementary sensitivity function
%   freq_Hz: frequency (in Hz) index of respT
%   C: transfer function of the feedback controller
%   Ts: sampling time
% 
% Xu Chen chx@uw.edu
% 2011-09-05
% 2015-02-08
% 2022-04-29
%
if nargin < 4
    Ts = 0.0004;
end
z = tf('z',Ts);
if nargin < 3
    C = 10000*(1+ 2*Ts/(1-z^(-1)) + 0.012* (1-z^(-1))/Ts);
end
if nargin < 2
    error('Not enough inputs')
end
fC = squeeze(freqresp(C,freq_Hz*2*pi));

fP = respT./(fC.*(1-respT));
gain = abs(fP);
phase = fixphasedata(angle(fP)*180/pi);
phase = mod(phase+180,360)-180; % edit phasewrap
function varargout = m_freq_resp_cal(y,u,Fs)
% function [Freq,mag,pha]=m_freq_resp_cal(y,u,Fs)
% This function computes the magnitude and phase of the frequency response between u (input)
% and y (output) using the spectral analysis.
% 
% X. Chen
% chx@uw.edu
% 2011-09-05
% 2022-04-29

FRF = [];

Ts = 1/Fs;
fmin_visu = 1;
fmax_visu = round(Fs/2);

Nfft = 2048/2;
window = hanning(Nfft);
noverlap = 0.75*Nfft;

[Txy,Freq] = tfestimate(u,y,window,noverlap,Nfft,Fs);% the unit of Freq is Hz

ind = find(Freq>fmin_visu & Freq<fmax_visu);

mag=20*log10(abs(Txy));
pha=180/pi*(phase(Txy));

figure;
subplot(211)
semilogx(Freq(ind),20*log10(abs(Txy(ind))));
if 0
    %%
    figure
    plot(Freq(ind),mag(ind))
    
    figure, bode(tf(1,[1 0] ))
end
grid
title('Frequency response')
xlim([min(Freq(ind)),max(Freq(ind))])
xl = xlim;
ylabel('Magnitude (dB)');
subplot(212)
semilogx(Freq(ind),pha(ind));grid
xlim(xl)
ylabel('Phase (degree)');
xlabel('Frequency (Hz)');

if nargout == 2
    varargout{1} = mag;
    varargout{2} = Freq;
elseif nargout == 3
    varargout{1} = mag;
    varargout{2} = Freq;
    varargout{3} = pha;
elseif nargout == 4
    varargout{1} = mag;
    varargout{2} = Freq;
    varargout{3} = pha;
    varargout{4} = Txy;
end
% This script provides an example of performing system identification using
% sinusoidal and random excitations
% associated file(s): "simuID.slx" simulink blocks
% 
% required functions: add "macs-matlab-toolbox-master" (available at https://github.com/macs-lab/macs-matlab-toolbox to the Matlab path;
%
% X. Chen
% chx@uw.edu
% 2015-02-4
% 2022-04-29

% clear
% clear all

% Sampling frequency
Fs = 1024; % change to the actual sampling frequency
% Generate a set of frequencies to perform the sinusoidal based system ID
initFreq = 1; % in Hz
endFreq = 100;
freqPoints = 20;
FreqRange = logspace(log10(initFreq),log10(endFreq),freqPoints);
% amplitude
amp_x = 2*ones(size(FreqRange))'; % change the amplitude to proper values in practice
% for data storage
tCell = cell(size(FreqRange));
xCell = cell(size(FreqRange));
yCell = xCell;
BODE_FREQ = zeros(freqPoints,1);
BODE_ARRAY = zeros(freqPoints,1);
%% system id based on sinusoidal excitation 
for ii = 1:freqPoints
    StopTime = 1/FreqRange(ii)*40; % 40 periods of the sinusoidal
    
    t = 0:1/Fs:StopTime;
    t = t(:);
    tCell{ii} = t;
    % number of points
    npts = length(t);
    
    % Create a sine wave 
    x = amp_x(ii) * (sin(2*pi*t*FreqRange(ii)));
    x = x(:);
    xCell{ii} = x;
    %% now feed the generated x file as the input to the system
    % the following data are dummies for simulation purpose
    % make proper CHANGEs to adapt to actual experiments
    num = [1];
    den = [1 0.5];
    sim('simuID')
    yCell{ii} = y;    
    % progressbar(min(1,f/(freqPoints+2)));
    
    %%
    % remove bias
    x = x - mean(x);
    y = y - mean(y);
    if 0
        %%
        % plot the signal
        figure(1)
        plot(t,x,t,y);
        xlabel('Time (s)');
        ylabel('Amplitude');
        legend('Signal x(t):sin(2*pi*t*20)',...
            'Signal y(t):sin(2*pi*t*20+lag)');
    end
    %%
%     tfestimate(x,y,[],[],[],Fs)
    % X,Y,WINDOW,NOVERLAP,NFFT,Fs
    [Txy,w] = tfestimate(x,y,[],[],[],Fs);
    BODE_FREQ(ii) = FreqRange(ii);
    try
        BODE_ARRAY(ii) = Txy(w==FreqRange(ii));
    catch
        if 0
            [ampTemp,phaTemp] = twoSine_amplitudeRatio_PhaseLag_direct(x,y,FreqRange(ii)*2*pi);
        else
            [ampTemp,phaTemp] = twoSine_amplitudeRatio_PhaseLag(x,y);
        end
        BODE_ARRAY(ii) = ampTemp * exp(1i*phaTemp);
%         Txy(w==FreqRange(ii))
    end
end

%% now do a random excitation
format long g
flag_spec = 0;
flag_fprint = 1;

ValUinit = 0;
ValAmpli = 1;
ValDecal = 0;
ValLgReg = 8;
ValDivi = 1;
Nsamp = 20*2^ValLgReg;
Tappli = 0;
%  "Entry parameters" are :
% 	  ValUinit  : Initial steady state
%     ValAmpli  : Magnitude
%     ValDecal  : Add-on DC component
%     ValLgReg  : Register length
%     ValDivi   : Frequency divider
%     Nsamp     : Number of samples
%     Tappli    : Application instant
prbs = create_prbs(ValUinit, ValAmpli, ValDecal, ValLgReg, ValDivi, Nsamp, Tappli);
prbs = prbs';
figure, plot(prbs)
if 0
    %%
    specPRBs = specCale(prbs,Fs);
    figure, plot(specPRBs.f,specPRBs.amp)
end
% different scaling of the random signal
prbs_8_1 = prbs;
prbs_8_1_scal1 = prbs*5e-7;
prbs_8_1_scal2 = prbs*10e-7;
prbs_8_1_scal50 = prbs*25e-7;

if 0
    %%
    save2tx(prbs_8_1_scal1,'prbs_8_1_scal1.txt');
    save2tx(prbs_8_1_scal2,'prbs_8_1_scal2.txt');
    save2tx(prbs_8_1_scal50,'prbs_8_1_scal50.txt');
end
if 0
    %%
    prbs_8_1_scal1 = [prbs_8_1_scal1,zeros(length(prbs_8_1_scal1),2)] ;
    prbs_8_1_scal2 = [prbs_8_1_scal2,zeros(length(prbs_8_1_scal2),2)] ;
    save('prbs_8_1_scal1.txt', 'prbs_8_1_scal1', '-ASCII','-DOUBLE')
    save('prbs_8_1_scal2.txt', 'prbs_8_1_scal2', '-ASCII')
end
% start with a small amplitude for safety, gradually incease the amplitude 
% in practice to find the best signal-to-noise ratio
x = prbs_8_1_scal1;
t = 0:1/Fs:(length(x)-1)/Fs;
t = t(:);
StopTime = t(end);
%%
sim('simuID')
y_prbs_8_1_scal1 = y;
%%
% figure, stairs(y_prbs_8_1_scal1(:,1))
% figure, stairs(prbs_8_1_scal1(:,1))
x = x - mean(x);
y = y - mean(y);

x = x(500:end);
y = y(500:end);
if 0
    %%
    figure,subplot(211),stairs(x),title('input')
    xlim([1,length(x)])
    subplot(212),stairs(y),title('output')
    xlim([1,length(y)])
end

if 0
    %%
    hgsave(gcf,'IOdata')
    saveas(gcf, 'IOdata', 'eps')
end
[mag,freq,pha] = m_freq_resp_cal(y,x,Fs);
if size(mag,1)<size(mag,2)
    mag = mag';
end
if size(freq,1)<size(freq,2)
    freq = freq';
end
% combine the sinusoidal sys id data with the random sys id data
BODE_FREQ = [FreqRange';freq(freq>max(FreqRange))];
BODE_MAG = [20*log10(abs(BODE_ARRAY(:)));mag(freq>max(FreqRange))];
BODE_PHA = [180/pi*(phase(BODE_ARRAY(:)));pha(freq>max(FreqRange))];
%%
figure, xbodeplot(tf(num,den,1/Fs))
% 
% figure;
subplot(211)
semilogx(BODE_FREQ,BODE_MAG,'r--','linewidth',2);
legend('actual','identified','location','best')
% grid
% title('Frequency response')
% xlim([min(BODE_FREQ),max(BODE_FREQ)])
% xl = xlim;
% ylabel('Magnitude (dB)');
subplot(212)
semilogx(BODE_FREQ,BODE_PHA,'r--','linewidth',2);
% grid
% xlim(xl)
% ylabel('Phase (degree)');
% xlabel('Frequency (Hz)');
if 0 % a second prbs profile that is longer in duration
    %%
    ValUinit = 0;
    ValAmpli = 1;
    ValDecal = 0;
    ValLgReg = 10;
    ValDivi = 1;
    Nsamp = 10*2^ValLgReg;
    Tappli = 0;
    %  "Entry parameters" are :
    % 	  ValUinit  : Initial steady state
    %     ValAmpli  : Magnitude
    %     ValDecal  : Add-on DC component
    %     ValLgReg  : Register length
    %     ValDivi   : Frequency divider
    %     Nsamp     : Number of samples
    %     Tappli    : Application instant
    prbs = create_prbs(ValUinit, ValAmpli, ValDecal, ValLgReg, ValDivi, Nsamp, Tappli);
    prbs = prbs';
    figure, plot(prbs)
    if flag_spec
        specPRBs = specCale(prbs,1);
        figure, plot(specPRBs.f,specPRBs.amp)
    end
    prbs_10_1 = prbs;
    prbs_10_1_scal1 = prbs*5e-7;
    prbs_10_1_scal2 = prbs*10e-7;
    prbs_10_1_scal50 = prbs*25e-7;
    
    if flag_fprint
        save2tx(prbs_10_1_scal1,'prbs_10_1_scal1.txt');
        save2tx(prbs_10_1_scal2,'prbs_10_1_scal2.txt');
        save2tx(prbs_10_1_scal50,'prbs_10_1_scal50.txt');
    else
        
        prbs_10_1_scal1 = [prbs_10_1_scal1,zeros(length(prbs_10_1_scal1),2)] ;
        prbs_10_1_scal2 = [prbs_10_1_scal2,zeros(length(prbs_10_1_scal2),2)] ;
        
        save('prbs_10_1_scal1.txt', 'prbs_10_1_scal1', '-ASCII')
        save('prbs_10_1_scal2.txt', 'prbs_10_1_scal2', '-ASCII')
    end
end
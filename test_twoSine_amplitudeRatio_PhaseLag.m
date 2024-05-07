% clear
% clear all

% set the lag
lag = pi/2;

% amplitude for sinusodial 1
amp_x = 2;

% amplitude for sinusodial 2
amp_y = 4;

% Sampling frequency
Fs = 1024;

% Time vector of 0.5 second
t = 0:1/Fs:0.5*(Fs-1)/Fs;

% number of points
npts = length(t);

% Create a sine wave of 20 Hz.
x = amp_x * (sin(2*pi*t*20) + ...
0.25*randn(1,npts)) + 5;

% Create a sine wave of 20 Hz
% with a phase of pi/2.
y = amp_y * (sin(2*pi*t*20 + lag)...
+ 0.25*randn(1,npts));

% basis vectors
sinb = sin(2*pi*t*20);
cosb = cos(2*pi*t*20);

% remove bias
x = x - mean(x);
y = y - mean(y);

% plot the signal
figure(1)
plot(t,x,t,y);
xlabel('Time (s)');
ylabel('Amplitude');
legend('Signal x(t):sin(2*pi*t*20)',...
'Signal y(t):sin(2*pi*t*20+lag)');

%%
tfestimate(x,y,[],[],[],Fs)
% X,Y,WINDOW,NOVERLAP,NFFT,Fs
[Txy,w] = tfestimate(x,y,[],[],[],Fs);
Txy(w==20)
%% direct computation
reX = sum(x.*sinb);
imX = sum(x.*cosb);
AX = sqrt(reX^2+imX^2);
phaX = atan(imX/reX);

reY = sum(y.*sinb);
imY = sum(y.*cosb);
AY = sqrt(reY^2+imY^2);
phaY = atan(imY/reY);

mag = AY/AX
pha = phaY-phaX
%% full fft way
% take the FFT
X=fft(x);
Y=fft(y);

% Calculate the numberof unique points
NumUniquePts = ceil((npts+1)/2);

figure(2)
subplot(211);
f = (0:NumUniquePts-1)*Fs/npts;
plot(f,abs(X(1:NumUniquePts)));
title('X(f) : Magnitude response');
ylabel('|X(f)|')
subplot(212)
plot(f,abs(Y(1:NumUniquePts)));
title('Y(f) : Magnitude response')
xlabel('Frequency (Hz)');
ylabel('|Y(f)|')

figure(3)
subplot(211)
plot(f,angle(X(1:NumUniquePts)));
title('X(f) : Phase response');
ylabel('Phase (rad)');
subplot(212)
plot(f,angle(Y(1:NumUniquePts)));
title('Y(f) : Phase response');
xlabel('Frequency (Hz)');
ylabel('Phase (rad)');

% Determine the max value and max point.
% This is where the sinusoidal
% is located. See Figure 2.
[mag_x idx_x] = max(abs(X));
[mag_y idx_y] = max(abs(Y));

% determine the phase difference
% at the maximum point.
px = angle(X(idx_x));
py = angle(Y(idx_y));
phase_lag = py - px

% determine the amplitude scaling
amplitude_ratio = mag_y/mag_x
%%

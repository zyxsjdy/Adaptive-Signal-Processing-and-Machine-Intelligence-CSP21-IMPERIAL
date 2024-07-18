%% Initialization
clc;
clear;
close all;

%% Setting Parameters
Fs = 1;  % normalized sampling frequency
n = 1024;  % number of samples
t = (0:n-1)/Fs;  % sampling time

S_imp = [2 zeros(1,n-1)];  % impulse
S_sin = sin(2*pi*0.1*t);  % sinusoid

%% DTFT of the ACF (definition 1)
% For impulse signal
ACF_imp = xcorr(S_imp,'biased');  % calculate the autocorrelation function
PSD_imp_1 = abs(fftshift(fft(ACF_imp)));  % calculate PSD

% For sinusoid signal
[ACF_sin,lags] = xcorr(S_sin,'biased');  % calculate the autocorrelation function
PSD_sin_1 = abs(fftshift(fft(ACF_sin)));  % calculate PSD

f_1 = (-n+1:n-1)/(2*n) * Fs;  % shifted frequency, zero-centered range

%% (average) Signal Power Over Frequencies (definition 2)
PSD_imp_2 = abs(fftshift(fft(S_imp))).^2 / n;
PSD_sin_2 = abs(fftshift(fft(S_sin))).^2 / n;

f_2 = (-n/2:n/2-1)/n * Fs;  % shifted frequency, zero-centered range

%% Result
figure
subplot(2,1,1)
plot(lags, ACF_imp,'Color','r','LineWidth',1.5)
grid on
title('Covariance Sequence (ACF) of the Impulse Signal')
xlabel('Lags')
ylabel('ACF')
xlim([min(lags) max(lags)])
subplot(2,1,2)
plot(lags, ACF_sin,'Color','b','LineWidth',0.1)
grid on
title('Covariance Sequence (ACF) of the Sinusoid Signal')
xlabel('Lags')
ylabel('ACF')
xlim([min(lags) max(lags)])

figure
plot(f_1,PSD_imp_1,'Color','b','LineWidth',1.5)
hold on
grid on
plot(f_2,PSD_imp_2,'Color','r','LineWidth',1.5,'LineStyle','--')
title('Periodogram of the Impulse Signal')
text(2,7,{'A simple plot','from 1 to 10'})
legend('Definition 1','Definition 2')
xlabel('Normalised Frequency (× \pi rad/sample)')
ylabel('PSD');

figure
plot(f_1,PSD_sin_1,'Color','b','LineWidth',1.5)
hold on
grid on
plot(f_2,PSD_sin_2,'Color','r','LineWidth',1.5,'LineStyle','--')
title('Periodogram of the Sinusoid Signal')
legend('Definition 1','Definition 2')
xlabel('Normalised Frequency (× \pi rad/sample)')
ylabel('PSD')

tilefigs([0 0.5 1 1])
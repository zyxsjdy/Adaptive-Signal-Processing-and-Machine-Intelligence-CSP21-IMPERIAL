%% Initialization
clc;
clear;
close all;

%% Setting Parameters
Fs = 1;  % normalized sampling frequency
n = 1000;  % number of samples
t = (0:n-1)/Fs;  % sampling time

noise = randn(1,n);  % WGN
S_nsin = sin(2*pi*0.1*t) + randn(1,n);  % noisy sinusoidal signal
f_noise = filter([1 1 1],2,noise);  % filtered WGN

f = (-n+1:n-1)/(2*n) * Fs;  % shifted frequency, zero-centered range

%% Using the Biased and Unbiased Estimators
% WGN
[ACF_n_b,lags] = xcorr(noise,'biased');
PSD_n_b = real(fftshift(fft(ifftshift(ACF_n_b))));
ACF_n_ub = xcorr(noise,'unbiased');
PSD_n_ub = real(fftshift(fft(ifftshift(ACF_n_ub))));

% noisy sinusoidal signal
ACF_s_b = xcorr(S_nsin,'biased');
PSD_s_b = real(fftshift(fft(ifftshift(ACF_s_b))));
ACF_s_ub = xcorr(S_nsin,'unbiased');
PSD_s_ub = real(fftshift(fft(ifftshift(ACF_s_ub))));

% filtered WGN
ACF_f_b = xcorr(f_noise,'biased');
PSD_f_b = real(fftshift(fft(ifftshift(ACF_f_b))));
ACF_f_ub = xcorr(f_noise,'unbiased');
PSD_f_ub = real(fftshift(fft(ifftshift(ACF_f_ub))));

%% Result
figure
subplot(3,2,1)
plot(lags,ACF_n_ub,'LineWidth',1.5)
hold on
grid on
plot(lags,ACF_n_b,'LineWidth',1.5)
title('ACF Estimate of WGN')
legend('Unbiased','Biased')
ylabel('Magnitude')
xlim([min(lags) max(lags)])
subplot(3,2,3)
plot(lags,ACF_s_ub,'LineWidth',1.5)
hold on
grid on
plot(lags,ACF_s_b,'LineWidth',1.5)
title('ACF Estimate of noisy sinusoidal signal')
legend('Unbiased','Biased')
ylabel('Magnitude')
xlim([min(lags) max(lags)])
subplot(3,2,5)
plot(lags,ACF_f_ub,'LineWidth',1.5)
hold on
grid on
plot(lags,ACF_f_b,'LineWidth',1.5)
title('ACF Estimate of filtered WGN')
legend('Unbiased','Biased')
xlabel('Lags')
ylabel('Magnitude')
xlim([min(lags) max(lags)])

subplot(3,2,2)
plot(f,PSD_n_ub,'LineWidth',1.5)
hold on
grid on
plot(f,PSD_n_b,'LineWidth',1.5)
legend('Unbiased','Biased')
title('Correlogram of WGN')
ylabel('Magnitude')
subplot(3,2,4)
plot(f,PSD_s_ub,'LineWidth',1.5)
hold on
grid on
plot(f,PSD_s_b,'LineWidth',1.5)
legend('Unbiased','Biased')
title('Correlogram of noisy sinusoidal signal')
ylabel('Magnitude')
subplot(3,2,6)
plot(f,PSD_f_ub,'LineWidth',1.5)
hold on
grid on
plot(f,PSD_f_b,'LineWidth',1.5)
legend('Unbiased','Biased')
title('Correlogram of filtered WGN')
xlabel('Normalised Frequency (× \pi rad/sample)')
ylabel('Magnitude')


tilefigs([0 0.4 0.7 1])
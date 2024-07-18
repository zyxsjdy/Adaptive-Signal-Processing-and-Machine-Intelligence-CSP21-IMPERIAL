%% Initialization
clc;
clear;
close all;

%% Setting Parameters
load('RRI-DATA.mat')

xRRI1 = detrend(xRRI1-mean(xRRI1));  % remove mean and trend
xRRI2 = detrend(xRRI2-mean(xRRI2));
xRRI3 = detrend(xRRI3-mean(xRRI3));

Fs = 4;  % sampling frequency, get from trial extraction
n_f = 1024;  % number of DFT points

%% Standard Periodogram
PSD_sta_1 = pow2db(periodogram(xRRI1,hamming(length(xRRI1)),n_f,Fs));
PSD_sta_2 = pow2db(periodogram(xRRI2,hamming(length(xRRI2)),n_f,Fs));
[PSD_sta_3,f] = periodogram(xRRI3,hamming(length(xRRI3)),n_f,Fs);
PSD_sta_3 = pow2db(PSD_sta_3);
n = length(PSD_sta_3);

%% AR Spectrum Estimate
p = 40;  % model orders
[ARcoef_est,error] = aryule(xRRI1,p);  % Yule-Walker Method
[h,~] = freqz(sqrt(error),ARcoef_est,n,Fs);
PSD_est_1 = pow2db(abs(h).^2);

p = 14;
[ARcoef_est,error] = aryule(xRRI2,p);
[h,~] = freqz(sqrt(error),ARcoef_est,n,Fs);
PSD_est_2 = pow2db(abs(h).^2);

p = 25;
[ARcoef_est,error] = aryule(xRRI3,p);
[h,f_a] = freqz(sqrt(error),ARcoef_est,n,Fs);
PSD_est_3 = pow2db(abs(h).^2);

%% Result
figure
subplot(3,1,1)
plot(f,PSD_sta_1,'LineWidth',1.5)
hold on
grid on
plot(f_a,PSD_est_1,'LineWidth',1.5)
title('PSD by Standard Periodogram and AR Spectrum Estimate for RRI Data1')
legend('Standard','Order 40')
xlabel('frequency (Hz)')
ylabel('Magnitude (dB)')
ylim([-100 0])

subplot(3,1,2)
plot(f,PSD_sta_2,'LineWidth',1.5)
hold on
grid on
plot(f_a,PSD_est_2,'LineWidth',1.5)
title('PSD by Standard Periodogram and AR Spectrum Estimate for RRI Data2')
legend('Standard','Order 14')
xlabel('frequency (Hz)')
ylabel('Magnitude (dB)')
ylim([-100 0])

subplot(3,1,3)
plot(f,PSD_sta_3,'LineWidth',1.5)
hold on
grid on
plot(f_a,PSD_est_3,'LineWidth',1.5)
title('PSD by Standard Periodogram and AR Spectrum Estimate for RRI Data3')
legend('Standard','Order 25')
xlabel('frequency (Hz)')
ylabel('Magnitude (dB)')
ylim([-100 0])

tilefigs([0 0.4 0.4 1])
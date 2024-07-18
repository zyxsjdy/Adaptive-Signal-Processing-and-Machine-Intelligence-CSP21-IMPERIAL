%% Initialization
clc;
clear;
close all;

%% Setting Parameters
load('RRI-DATA.mat')

xRRI1 = detrend(xRRI1-mean(xRRI1));  % remove mean and trend
xRRI2 = detrend(xRRI2-mean(xRRI2));
xRRI3 = detrend(xRRI3-mean(xRRI3));

Fs = 4;  % sampling frequency
n_f = 1024;  % number of DFT points

%% Standard Periodogram
PSD_sta_1 = pow2db(periodogram(xRRI1,hamming(length(xRRI1)),n_f,Fs));  % apply the satandard periodogram
PSD_sta_2 = pow2db(periodogram(xRRI2,hamming(length(xRRI2)),n_f,Fs));
PSD_sta_3 = pow2db(periodogram(xRRI3,hamming(length(xRRI3)),n_f,Fs));

%% Averaged Periodogram
PSD_ave_11 = pow2db(pwelch(xRRI1,hamming(50*Fs),0,n_f,Fs));  % apply the averaged periodogram
PSD_ave_12 = pow2db(pwelch(xRRI1,hamming(150*Fs),0,n_f,Fs));

PSD_ave_21 = pow2db(pwelch(xRRI2,hamming(50*Fs),0,n_f,Fs));
PSD_ave_22 = pow2db(pwelch(xRRI2,hamming(150*Fs),0,n_f,Fs));

PSD_ave_31 = pow2db(pwelch(xRRI3,hamming(50*Fs),0,n_f,Fs));
[PSD_ave_32,f] = pwelch(xRRI3,hamming(150*Fs),0,n_f,Fs);  % f is the corresponding normalized frequency
PSD_ave_32 = pow2db(PSD_ave_32);

%% Result
figure
subplot(3,1,1)
plot(f,PSD_sta_1,'LineWidth',1.5)
hold on
grid on
plot(f,PSD_ave_11, 'LineWidth',1.5)
plot(f,PSD_ave_12, 'LineWidth',1.5)
title('PSD of the RRI Data1')
legend('Standard','\Deltat = 50','\Deltat = 150')
xlabel('frequency (Hz)')
ylabel('Magnitude (dB)')
ylim([-100 0])

subplot(3,1,2)
plot(f,PSD_sta_2,'LineWidth',1.5)
hold on
grid on
plot(f,PSD_ave_21, 'LineWidth',1.5)
plot(f,PSD_ave_22, 'LineWidth',1.5)
title('PSD of the RRI Data2')
legend('Standard','\Deltat = 50','\Deltat = 150')
xlabel('frequency (Hz)')
ylabel('Magnitude (dB)')
ylim([-100 0])

subplot(3,1,3)
plot(f,PSD_sta_3,'LineWidth',1.5)
hold on
grid on
plot(f,PSD_ave_31, 'LineWidth',1.5)
plot(f,PSD_ave_32, 'LineWidth',1.5)
title('PSD of the RRI Data3')
legend('Standard','\Deltat = 50','\Deltat = 150')
xlabel('frequency (Hz)')
ylabel('Magnitude (dB)')
ylim([-100 0])

tilefigs([0 0.4 0.5 1])
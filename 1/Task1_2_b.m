%% Initialization
clc;
clear;
close all;

%% Setting Parameters
load('EEG_Data_Assignment1.mat');
Fs = fs;  % sampling frequency
s = POz;  % signal
s1 = detrend(s - mean(s));  % data remove mean and trend

n = length(s1);  % number of samples
n_f = 5 * Fs;  % number of DFT points, 5 per Hz
L = [10 5 1];  % window lengths

%% Standard Periodogram Approach
[PSD_sta,f] = periodogram(s1,hamming(n),n_f,Fs);
PSD_sta = pow2db(PSD_sta);

%% Averaged Periodogram Approach
% Welch's power spectral density estimate using hamming window
PSD_ave_1 = pow2db(pwelch(s1,hamming(L(1)*Fs),0,n_f,Fs));
PSD_ave_2 = pow2db(pwelch(s1,hamming(L(2)*Fs),0,n_f,Fs));
PSD_ave_3 = pow2db(pwelch(s1,hamming(L(3)*Fs),0,n_f,Fs));

%% Result
figure
plot(f,PSD_sta,'LineWidth',1.5)
grid on
hold on
plot(f,PSD_ave_1,'LineWidth',1.5)
title('Standard Periodogram Approach')
legend('Standard','\Deltat = 10s')
xlabel('frequency (Hz)')
ylabel('Magnitude (dB)')
xlim([0 60])

figure
plot(f,PSD_ave_1,'LineWidth',1.5,'Color','r')
hold on
plot(f,PSD_ave_2,'LineWidth',1.5,'Color','b')
plot(f,PSD_ave_3,'LineWidth',1.5,'Color','g')
grid on
title('Averaged Periodogram Approach')
legend('\Deltat = 10s','\Deltat = 5s','\Deltat = 1s')
xlabel('frequency (Hz)')
ylabel('Magnitude (dB)')
xlim([0 60])

figure
plot(f,PSD_sta,'LineWidth',1.5)
hold on
plot(f,PSD_ave_1,'LineWidth',1.5)
plot(f,PSD_ave_2,'LineWidth',1.5)
plot(f,PSD_ave_3,'LineWidth',1.5)
grid on
legend('Standard','\Deltat = 10s','\Deltat = 5s','\Deltat = 1s')
xlabel('frequency (Hz)')
ylabel('Magnitude (dB)')
xlim([0 60])

tilefigs([0 0.5 1 1])
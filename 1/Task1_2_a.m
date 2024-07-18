%% Initialization
clc;
clear;
close all;

%% Setting Parameters
load sunspot.dat
Fs = 1;  % normalized sampling frequency
n = length(sunspot);  % number of samples
t = sunspot(:,1);  % sampling time

% Preprocessing
s = sunspot(:,2) + 1e-1;  % sunspot time series

%% Periodogram-based Spectral Estimation
% Welch's power spectral density estimate using hamming window
[PSD1,f] = pwelch(s,hamming(n),0,n,Fs);

% PSD (remove mean and trend)
s_r = detrend(s-mean(s));
PSD2 = pwelch(s_r,hamming(n),0,n,Fs);

% PSD (log and remove mean)
s_l = log(s) - mean(log(s));
PSD3 = pwelch(s_l,hamming(n),0,n,Fs);

%% Result
figure
plot(t,s,'LineWidth',1.5)
hold on
grid on
plot(t,s_r,'LineWidth',1.5)
plot(t,s_l,'LineWidth',1.5)
title('Sunspot Time Series')
legend('Original','Remove mean & trend','Log & Remove mean')
xlabel('Year')
ylabel('Sunspots')
xlim([min(t) max(t)])

figure
sgtitle('Periodogram-based Spectral Estimation for Sunspot Time Series')
subplot(1,3,1)
plot(f,PSD1,'LineWidth',1.5)
hold on
grid on
plot(f,PSD2,'LineWidth',1.5,'LineStyle','--')
legend('Original','Remove Mean & Trend')
xlabel('frequency (cycles/year)')
ylabel('PSD')
xlim([min(f) max(f)])

subplot(1,3,2)
plot(f,PSD2,'LineWidth',1.5)
grid on
legend('Remove Mean & Trend')
xlabel('frequency (cycles/year)')
ylabel('PSD')
axis([min(f) max(f) 0 1.1*max(PSD2)])

subplot(1,3,3)
plot(f,PSD3,'LineWidth',1.5)
grid on
legend('Log & Remove Mean')
xlabel('frequency (cycles/year)')
ylabel('PSD')
axis([min(f) max(f) 0 1.1*max(PSD3)])


tilefigs([0 0.6 0.7 1])
%% Initialization
clc;
clear;
close all;

%% Setting Parameters
Fs = 1;  % normalized sampling frequency
n = 1000;  % number of samples
t = (0:n-1)/Fs;  % sampling time

S = exp(1i*2*pi*0.3*t) + exp(1i*2*pi*0.32*t);  % generate the signal with two complex exponentials
L = 30;  % Signal length
S = S(1:L);

%% PSD of Signals of Different Length
loop = 100;
PSD = zeros(256,loop);
for j = 1:loop
    noise = sqrt(0.2/2) * (randn(1,L) + 1i*randn(1,L));  % generate noise
    S_n = S + noise;

    [X,R] = corrmtx(S_n,14,'mod');  % Return an autocorrelation estimate for the input signal
    [PSD(:,j),F] = pmusic(R,2,[],Fs,'corr');  % Estimate spectrum using MUSIC algorithm
end

PSD_mean = mean(PSD,2);  % calculate the mean
PSD_std = std(PSD,[],2);  % calculate the standard deviation

%% Result
figure
subplot(1,2,1)
plot(F,PSD(:,1),'Color','c','LineWidth',1.5)
set(gca,'xlim',[0.25 0.40])
hold on
grid on
for j = 2:loop
    plot(F,PSD(:,j),'Color','c','LineWidth',1.5)
end
plot(F,PSD_mean,'Color','b','LineWidth',1)
title('PSD using MUSIC Algorithm')
xlabel('Normalised Frequency (× \pi rad/sample)')
ylabel('Magnitude')

subplot(1,2,2)
plot(F,PSD_std,'LineWidth',1.5)
set(gca,'xlim',[0.25 0.40])
grid on
title('Standard Deviation of PSD using MUSIC Algorithm')
xlabel('Normalised Frequency (× \pi rad/sample)')
ylabel('Magnitude')

tilefigs([0 0.4 0.7 1])
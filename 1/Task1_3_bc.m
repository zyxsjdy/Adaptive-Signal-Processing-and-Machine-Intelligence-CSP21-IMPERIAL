%% Initialization
clc;
clear;
close all;

%% Setting Parameters
Fs = 1;  % normalized sampling frequency
n = 100;  % number of samples
t = (0:n-1)/Fs;  % sampling time

S = sin(2*pi*0.1*t);  % sinusoidal signal

f = (-n+1:n-1)/(2*n) * Fs;  % shifted frequency, zero-centered range

%% Generate the PSD Estimate of Several Realisations of a Random Process
loop = 100;
PSD = zeros(loop,2*n-1);
for i = 1:loop
    S_n = S + 0.2*randn(1,n);  % noisy sinusoidal signal
    [ACF,~] = xcorr(S_n,'biased');
    PSD(i,:) = real(fftshift(fft(ifftshift(ACF))));
end

PSD_mean = mean(PSD);  % calculate the mean
PSD_std = std(PSD);  % calculate the standard deviation

% dB version
PSD_dB = pow2db(PSD);
PSD_mean_dB = mean(PSD_dB);
PSD_std_dB = std(PSD_dB);

%% Result
figure
subplot(1,2,1)
plot(f,PSD(1,:),'Color','c','LineWidth',1.5)
hold on
for i = 2:loop
    plot(f,PSD(i,:),'Color','c','LineWidth',1.5)
end
grid on
plot(f,PSD_mean,'Color','b','LineWidth',1)
title('PSD of Sinusoids Corrupted by WGN')
xlabel('Normalised Frequency (× \pi rad/sample)')
ylabel('Magnitude')
subplot(1,2,2)
plot(f,PSD_std,'LineWidth',1.5)
grid on
title('Standard Deviation of PSD of Sinusoids Corrupted by WGN')
xlabel('Normalised Frequency (× \pi rad/sample)')
ylabel('Magnitude')

figure
subplot(1,2,1)
plot(f,PSD_dB(1,:),'Color','c','LineWidth',1.5)
hold on
for i = 2:loop
    plot(f,PSD_dB(i,:),'Color','c','LineWidth',1.5)
end
grid on
plot(f,PSD_mean_dB,'Color','b','LineWidth',1)
title('PSD of Sinusoids Corrupted by WGN')
xlabel('Normalised Frequency (× \pi rad/sample)')
ylabel('Magnitude (dB)')
subplot(1,2,2)
plot(f,PSD_std_dB,'LineWidth',1.5)
grid on
title('Standard Deviation of PSD of Sinusoids Corrupted by WGN')
xlabel('Normalised Frequency (× \pi rad/sample)')
ylabel('Magnitude (dB)')

tilefigs([0 0.4 0.7 1])
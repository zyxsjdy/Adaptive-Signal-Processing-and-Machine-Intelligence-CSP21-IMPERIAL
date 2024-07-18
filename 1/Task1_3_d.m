%% Initialization
clc;
clear;
close all;

%% Setting Parameters
Fs = 1;  % normalized sampling frequency
n = 1000;  % number of samples
t = (0:n-1)/Fs;  % sampling time

f = (-n/2:n/2-1)/n * Fs;  % symmetrical frequency points

S = exp(1i*2*pi*0.3*t) + exp(1i*2*pi*0.32*t);  % generate the signal with two complex exponentials
noise = sqrt(0.2/2) * (randn(1,n) + 1i*randn(1,n));  % generate noise
S_n = S + noise;

L = [20 30 40 50 60];  % Signal length

%% PSD of Signals of Different Length
PSD = zeros(length(L),n);
for i = 1:length(L)
    s = S_n(1:L(i));  % take length L(i)
    PSD(i,:) = pow2db(abs(fftshift(fft(s,n))).^2 / L(i));  % calculate PSD
end

%% Result
figure
for i = 1:length(L)
    plot(f,PSD(i,:),'LineWidth',1.5)
    hold on
end
grid on
title('Periodogram of Signals of Different Length')
legend(int2str(L'))
xlabel('Normalised Frequency (× \pi rad/sample)')
ylabel('PSD (dB)')
axis([0.25 0.4 -30 25])

tilefigs([0 0.4 0.4 1])
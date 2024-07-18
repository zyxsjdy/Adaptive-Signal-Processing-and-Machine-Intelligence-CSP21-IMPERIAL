%% Initialization
clc;
clear;
close all;

%% Setting Parameters
n = 1500;  % number of samples
Fs = 1500;  % sampling frequency
L = 3;  % phase is composed of 3 pieces
n_L = 500;  % length of each piece
L_AR = 1;  % AR order

%% Find the AR(1) Coefficient
% Generate the frequency function
f = 100*ones(1,1500);
for i = 501:1500
    if i <= 1000
        f(i) = f(i) + (i-500)/2;
    else
        f(i) = f(i) + ((i-1000)/25)^2;
    end
end
phi = cumsum(f);  % find the phase
y = exp(1i*2*pi/Fs*phi) + sqrt(0.05/2)*(randn(1,n)+1i*randn(1,n));  % FM signal

AR_coef = aryule(y,L_AR);  % using the Yule-Walker method to estimate the AR parameter
[h,w] = freqz(1,AR_coef,n,Fs);
PSD = pow2db(abs(h).^2);  % PSD of the signal

% Estimation for each pieces
w3 = zeros(n_L,L);
PSD3 = zeros(n_L,L);
for i = 1:L
    AR_coef = aryule(y(1+n_L*(i-1):n_L*i),L_AR);
    [h,w3(:,i)] = freqz(1,AR_coef,n_L,Fs);
    PSD3(:,i) = pow2db(abs(h).^2);
end

%% Result
figure
subplot(2,1,1)
plot(f)
grid on
title('Frequency of the FM Signal')
xlabel('Sample')
ylabel('frequency (Hz)')
subplot(2,1,2)
plot(angle(exp(1i*2*pi/Fs*phi)))
grid on
title('Phase of the FM signal')
xlabel('Sample')
ylabel('frequency (Hz)')

figure
sgtitle('PSD of the FM Signal using AR(1) Model')

subplot(L+1,1,1)
plot(w,PSD,'LineWidth',1.5)
grid on
xlabel('frequency (Hz)')
ylabel('Magnitude (dB)')
for i = 1:L
    subplot(L+1,1,i+1)
    plot(w3(:,i),PSD3(:,i),'LineWidth',1.5)
    grid on
    title(['Segment' num2str(i)])
    xlabel('frequency (Hz)')
    ylabel('Magnitude (dB)')
end

tilefigs([0 0.4 0.7 1])
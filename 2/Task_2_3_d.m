%% Initialization
clc;
clear;
close all;

%% Setting Parameters
load('EEG_Data_Assignment2.mat');
POz = POz';

n = length(POz);  % number of samples
t = (0:n-1)/fs;  % sampling time
n_window = 5*fs;  % window length
n_overlap = round(0.5*n_window);  % window overlap length
n_FFT = 15*fs;  % FFT points

u = [0.1,0.01,0.001];
n_u = length(u);
M = [3 10 15];
n_M = length(M);

%% Adaptive noise cancellation
S_nsin = sin(2*pi*50*t) + sqrt(0.05) * randn(1,n);
POz_delay = [0 POz(1:end-1)];

S_ANC = zeros(1,n,n_M,n_u);
MSPE_ANC = zeros(n_M,n_u);
for i = 1:n_M
    S_delay = zeros(M(i),n);
    for k = 1:M(i)
        S_delay(k,:) = [zeros(1,k) S_nsin(1:n-k)];  % delay the noisy signal for ANC
    end
    for j = 1:n_u
        [~,noise] = LMS_AR_1(S_delay,POz_delay,u(j));  % apply LMS algorithm to estimate the noise
        S_ANC(:,:,i,j) = POz_delay - noise;  % substract the estimated noise to get the ANC signal
        MSPE_ANC(i,j) = sum((POz-S_ANC(:,:,i,j)).^2)/n;  % acculumate the squared error
    end
end
PSD_POz = periodogram(POz,hamming(n),n_FFT,fs);
[PSD_ANC,f] = periodogram(S_ANC(:,:,2,2),hamming(n),n_FFT,fs);

%% Result
figure
plot(M,MSPE_ANC(:,1),'LineWidth',1.5)
hold on
grid on
plot(M,MSPE_ANC(:,2),'LineWidth',1.5)
plot(M,MSPE_ANC(:,3),'LineWidth',1.5)
title('MSPE vs Filter Order')
legend('\mu = 0.1','\mu = 0.01','\mu = 0.001')
xlabel('Filter Order')
ylabel('Magnitude')
xlim([min(M) max(M)])

figure
subplot(2,1,1)
plot(f,pow2db(PSD_POz),'LineWidth',1.5)
hold on
grid on
plot(f,pow2db(PSD_ANC))
title('Periodogram of Original POz and ANC Result')
legend('POz','ANC')
ylabel('Magnitude (dB)')
ylim([-160 -90])
subplot(2,1,2)
plot(f,pow2db(abs((PSD_POz-PSD_ANC).^2)),'LineWidth',1.5)
grid on
title('Difference between Original POz and ANC Result')
xlabel('frequency (Hz)')
ylabel('Magnitude (dB)')


figure
spectrogram(POz,n_window,n_overlap,n_FFT,fs,'yaxis')
title('Spectrogram of POz')
ylim([0 60])

for i = 1:n_M
    figure
    sgtitle(['Spectrogram of Denoised POz, M = ' num2str(M(i))])
    for j = 1:n_u
        subplot(1,n_u,j)
        spectrogram(S_ANC(:,:,i,j),n_window,n_overlap,n_FFT,fs,'yaxis')
        title(['\mu = ' num2str(u(j))])
        ylim([0 60])
    end
end

tilefigs([0 0.05 0.95 1])
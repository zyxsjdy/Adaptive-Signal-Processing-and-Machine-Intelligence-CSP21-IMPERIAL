%% Initialization
clc;
clear;
close all;

%% Setting Parameters
n = 1000;  % number of samples
t = 0:n-1;  % sampling time
S_sin = sin(2*pi*0.01*t);  % sinusoidal signal
loop = 100;
n_var = 1;  % innovation variance
MA_coef = [0 0.5];  % MA coefficient
u = 0.01;  % step size
M = 2;  % filter length
delay = 3;  % delay

%% Adaptive line enhancer and adaptive noise cancellation
model = arima('Constant',0,'MA',MA_coef,'Variance',n_var);  % Create univariate autoregressive integrated moving average (ARIMA) model
[S_MA,~,~] = simulate(model,n,'NumPaths',loop);  % Simulate sample paths of the model
S_MA = S_MA';  % MA Signal, considered as the noise

noise_2nd = 0.7*S_MA + 0.05*randn(size(S_MA));  % secondary noise
corr_ = mean(diag(corr(S_MA',noise_2nd')));  % average correlation coefficient between S_MA and noise_2nd
fprintf(['correlation coefficient between S_MA and noise_2nd: ' num2str(corr_) '\n'])

S_nsin = zeros(loop,n);
N_delay = zeros(M,n);
S_ALE = zeros(loop,n);
S_ANC = zeros(loop,n);
e_ALE = 0;
e_ANC = 0;
e_ANC_500 = 0;
for i = 1:loop
    S_nsin(i,:) = S_sin + S_MA(i,:);  % mix the clean signal and the noise
    [~,S_ALE(i,:),~] = LMS_AR(S_nsin(i,:),M,delay,u,0);  % apply LMS algorithm
    e_ALE = e_ALE + sum((S_sin-S_ALE(i,:)).^2)/n;  % acculumate the squared error

    S_delay = [0 S_nsin(i,1:end-1)];  % delay the signal for ANC
    for j = 1:M
        N_delay(j,:) = [zeros(1,j) noise_2nd(i,1:n-j)];  % delay the 2nd noise for ANC
    end
    [~,noise] = LMS_AR_1(N_delay,S_delay,u);  % apply LMS algorithm to estimate the noise
    S_ANC(i,:) = S_delay - noise;  % substract the estimated noise to get the ANC signal
    e_ANC = e_ANC + sum((S_sin-S_ANC(i,:)).^2)/n;  % acculumate the squared error
    e_ANC_500 = e_ANC_500 + sum((S_sin(501:end)-S_ANC(i,501:end)).^2)/n;  % acculumate the squared steady state error
end
S_ALE_ave = mean(S_ALE,1);
S_ANC_ave = mean(S_ANC,1);
MSPE_ALE = e_ALE/loop;
MSPE_ANC = e_ANC/loop;
MSPE_ANC_500 = e_ANC_500/loop;

fprintf(['MSPE_ALE: ' num2str(MSPE_ALE) '\n'])
fprintf(['MSPE_ANC: ' num2str(MSPE_ANC) '\n'])
fprintf(['MSPE_ANC in the last 500 samples: ' num2str(MSPE_ANC_500) '\n'])

%% Result
figure
subplot(3,1,1)
plot(t,S_nsin(1,:),'Color','b','LineWidth',1.5)
hold on
grid on
plot(t,S_ALE(1,:),'Color','r','LineWidth',1.5)
plot(t,S_sin,'Color','k','LineWidth',1.5)
for i = 1:loop
    plot(t,S_nsin(i,:),'Color','b','LineWidth',1.5)
end
for i = 1:loop
    plot(t,S_ALE(i,:),'Color','r','LineWidth',1.5)
end
plot(t,S_sin,'Color','k','LineWidth',1.5)
title(['ALE using LMS Algorithm, \Delta = 3, M = 2, MSPE = ' num2str(MSPE_ALE)])
legend('Noise-corrpupted Signal','Estimated De-noised Signal','Original Signal')
xlabel('Sample')
ylabel('Amplitude')

subplot(3,1,2)
plot(t,S_nsin(i,:),'Color','b','LineWidth',1.5)
hold on
grid on
plot(t,S_ANC(i,:),'Color','r','LineWidth',1.5)
plot(t,S_sin,'Color','k','LineWidth',1.5)
for i = 1:loop
    plot(t,S_nsin(i,:),'Color','b','LineWidth',1.5)
end
for i = 1: loop
    plot(t,S_ANC(i,:),'Color','r','LineWidth',1.5)
end
plot(t,S_sin,'Color','k','LineWidth',1.5)
title(['ANC using LMS Algorithm, \Delta = 3, M = 2, MSPE = ' num2str(MSPE_ANC)])
legend('Noise-corrpupted Signal','Estimated De-noised Signal','Original Signal')
xlabel('Sample')
ylabel('Amplitude')

subplot(3,1,3)
plot(t,S_sin,'Color','k','LineWidth',1.5)
hold on
grid on
plot(t,S_ALE_ave,'Color','b','LineWidth',1.5)
plot(t,S_ANC_ave,'Color','r','LineWidth',1.5)
title('Averaged 100 Realizations')
legend('Original Signal','ALE Result','ANC Result')
xlabel('Sample')
ylabel('Amplitude')

tilefigs([0 0.4 0.4 1])
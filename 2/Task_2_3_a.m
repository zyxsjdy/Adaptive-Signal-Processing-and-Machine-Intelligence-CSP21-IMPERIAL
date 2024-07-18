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
M = 5;  % filter length
delay = 7;  % maximum delay

%% Adaptive Line Enhancer (ALE)
model = arima('Constant',0,'MA',MA_coef,'Variance',n_var);  % Create univariate autoregressive integrated moving average (ARIMA) model
[S_MA,~,~] = simulate(model,n,'NumPaths',loop);  % Simulate sample paths of the model
S_MA = S_MA';  % MA Signal, considered as the noise

S_nsin = zeros(1,n,loop,delay);
S_ALE = zeros(1,n,loop,delay);
MSPE = zeros(delay,1);
for j = 1:delay
    error = 0;
    for i = 1:loop
        S_nsin(:,:,i,j) = S_sin + S_MA(i,:);  % mix the clean signal and the noise
        [~,S_ALE(:,:,i,j),~] = LMS_AR(S_nsin(:,:,i,j),M,j,u,0);  % apply LMS algorithm
        error = error + sum((S_sin-S_ALE(:,:,i,j)).^2)/n;  % acculumate the squared error
    end
    MSPE(j) = error/loop;  % mean
end

%% Result
figure
sgtitle('ALE using LMS Algorithm, filter length M = 5')

for j = 1:delay
    subplot(1,delay,j)
    plot(t,S_nsin(:,:,1,j),'Color','b','LineWidth',1.5)
    hold on
    grid on
    plot(t,S_ALE(:,:,1,j),'Color','r','LineWidth',1.5)
    plot(t,S_sin,'Color','k','LineWidth',1.5)
    for i = 1:loop
        plot(t,S_nsin(:,:,i,j),'Color','b','LineWidth',1.5)
    end
    for i = 1:loop
        plot(t,S_ALE(:,:,i,j),'Color','r','LineWidth',1.5)
    end
    plot(t,S_sin,'Color','k','LineWidth',1.5)
    title(['\Delta = ' num2str(j)])
    xlabel('Sample')
    ylabel('Magnitude')
end
legend('Noise-corrpupted Signal','Estimated De-noised Signal','Original Signal')

figure
plot(MSPE,'LineWidth',1.5)
grid on
title('Mean Square Prediction Error')
xlabel('\Delta')
ylabel('Magnitude')

tilefigs([0 0.4 0.7 1])
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
M = 1:25;  % filter length
n_M = length(M);
delay = 1:25;
n_delay = length(delay);

%% Adaptive Line Enhancer (ALE)
model = arima('Constant',0,'MA',MA_coef,'Variance',n_var);  % Create univariate autoregressive integrated moving average (ARIMA) model
[S_MA,~,~] = simulate(model,n,'NumPaths',loop);  % Simulate sample paths of the model
S_MA = S_MA';  % MA Signal, considered as the noise

MSPE = zeros(n_M,n_delay);
for k = 1:n_M
    for j = 1:n_delay
        error = 0;
        for i = 1:loop
            S_nsin = S_sin + S_MA(i,:);  % mix the clean signal and the noise
            [~,S_ALE,~] = LMS_AR(S_nsin,M(k),j,u,0);  % apply LMS algorithm
            error = error + sum((S_sin-S_ALE).^2)/n;  % acculumate the squared error
        end
        MSPE(k,j) = error/loop;  % mean
    end
end

%% Result
figure
subplot(1,2,1)
plot(delay,MSPE(5,:),'LineWidth',1.5)
hold on
grid on
plot(delay,MSPE(10,:),'LineWidth',1.5)
plot(delay,MSPE(15,:),'LineWidth',1.5)
plot(delay,MSPE(20,:),'LineWidth',1.5)
title('MSPE vs Delay')
legend('M = 5','M = 10','M = 15','M = 20')
xlabel('Delay')
ylabel('Magnitude')
xlim([3 max(delay)])

subplot(1,2,2)
plot(M,MSPE(:,3),'LineWidth',1.5)
grid on
title('MSPE vs Filter Order, \Delta = 3')
xlabel('Filter Length')
ylabel('Magnitude')
xlim([1 max(M)])

tilefigs([0 0.4 0.7 1])
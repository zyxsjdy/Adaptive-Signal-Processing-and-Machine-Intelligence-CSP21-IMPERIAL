%% Initialization
clc;
clear;
close all;

%% Setting Parameters
loop = 100;
n = 1000;  % number of samples
n_var = 0.25;  % noise variance
AR_coef = [0.1 0.8];  % AR coefficients
L_AR = length(AR_coef);
u = [0.05 0.01];  % step sizes
L_u = length(u);

%% LMS Adaptive Predictor
model = arima('Constant',0,'AR',AR_coef,'Variance',n_var);  % Create univariate autoregressive integrated moving average (ARIMA) model
S_AR = (simulate(model,n,'NumPaths',loop))';  % Simulate sample paths of the model

error_sqr = zeros(L_u,n,loop);
error_ave = zeros(L_u,n);
for k = 1:L_u
    for i = 1:loop
        S = S_AR(i,:);
        [~,~,error] = LMS_AR(S,L_AR,1,u(k),0);  % find the error
        error_sqr(k,:,i) = error.^2;
    end
    error_ave(k,:) = pow2db(mean(error_sqr(k,:,:),3));  % find the average error
end

%% Result
figure
subplot(2,1,1)
plot(pow2db(error_sqr(1,:,1).^2),'LineWidth',1.5)  % show the error of the first realisation as example
hold on
grid on
plot(pow2db(error_sqr(2,:,1).^2),'LineWidth',1.5)
title('Squared Prediction Error')
legend('\mu = 0.05','\mu = 0.01')
xlabel('Sample')
ylabel('Magnitude (dB)')

subplot(2,1,2)
plot(error_ave(1,:),'LineWidth',1.5)
hold on
grid on
plot(error_ave(2,:),'LineWidth',1.5)
title('Average Squared Prediction Error of 100 Realizations')
legend('\mu = 0.05','\mu = 0.01')
xlabel('Sample')
ylabel('Magnitude (dB)')

tilefigs([0 0.4 0.4 1])
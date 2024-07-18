%% Initialization
clc;
clear;
close all;

%% Setting Parameters
load('time-series.mat');
S = (y-mean(y))';  % remove the mean
p = 4;  % AR order
delay = 1;
u = 1e-5;  % learning rate

%% LMS
[~,S_est,error] = LMS_AR(S,p,delay,u,0);  % apply the LMS prediction
MSE = pow2db(mean(abs(error).^2));  % MSE
MSE_s = pow2db(mean(abs(error(501:end)).^2));  % MSE in steady state
R_p = pow2db(var(S_est)/var(error));  % prediction gain
R_p_s = pow2db(var(S_est(501:end))/var(error(501:end)));  % prediction gain in steady state

%% Result
fprintf(['Overall MSE = ' num2str(MSE) ' dB\n'])
fprintf(['Overall Prediction gain = ' num2str(R_p) ' dB\n'])
fprintf(['MSE in steady state = ' num2str(MSE_s) ' dB\n'])
fprintf(['Prediction gain in steady state = ' num2str(R_p_s) ' dB\n'])

figure
plot(S,'Linewidth',1.5)
hold on
grid on
plot(S_est,'Linewidth',1.5)
title('Zero-mean Version of y and its One-step Ahead Prediction using the LMS Algorithm')
legend('Zero-mean Version','One-step Ahead Prediction')
xlabel('Sample')
ylabel('Magnitude')

tilefigs([0 0.5 0.6 1])
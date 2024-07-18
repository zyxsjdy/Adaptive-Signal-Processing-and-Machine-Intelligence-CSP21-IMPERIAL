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

MSE = zeros(L_u,loop);
for k = 1:L_u
    for i = 1:loop
        S = S_AR(i,:);
        [~,~,error] = LMS_AR(S,L_AR,1,u(k),0);  % find the error
        MSE(k,i) = mean(error(501:end).^2);  % (501:end) remove the transient output
    end
end

EMSE = (mean(MSE-n_var,2))';  % Excess Mean Square Error
M_app = EMSE/n_var;  % misadjustment

Rxx = [25/27  25/54;
       25/54  25/27];  % autocorrelation matrix
M_the = u/2*trace(Rxx);  % approximated misadjustment

%% Result
fprintf(['misadjustment: ' num2str(M_app,4) '\n'])
fprintf(['approximated misadjustment: ' num2str(M_the,4) '\n'])

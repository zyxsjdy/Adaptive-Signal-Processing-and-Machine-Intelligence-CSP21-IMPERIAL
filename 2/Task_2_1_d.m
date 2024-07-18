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
w = zeros(L_AR,n,loop);  % weight
w_ave = zeros(L_AR,n,L_u);  % average weight
for k = 1:L_u
    for i = 1:loop
        S = S_AR(i,:);
        [w(:,:,i),~,~] = LMS_AR(S,L_AR,1,u(k),0);  % find the error
    end
    w_ave(:,:,k) = mean(w(:,:,:),3);
end

%% Result
u5 = [mean(w_ave(1,800:end,1)) mean(w_ave(2,800:end,1))];  % Coefficients for u = 0.05
u1 = [mean(w_ave(1,800:end,2)) mean(w_ave(2,800:end,2))];  % Coefficients for u = 0.01
fprintf(['coefficients with u = 0.05:     ' num2str(u5,4) '\n'])
fprintf(['coefficients with u = 0.01:     ' num2str(u1,4) '\n'])

figure
subplot(L_u,1,1)
plot(w_ave(1,:,1),'LineWidth',1.5)
hold on
grid on
plot(w_ave(2,:,1),'LineWidth',1.5)
plot([0 n],[AR_coef(1) AR_coef(1)],'LineWidth',1.5,'LineStyle','--')
plot([0 n],[AR_coef(2) AR_coef(2)],'LineWidth',1.5,'LineStyle','--')
title('Updating Curve of the Adaptive Filter Coefficients for \mu = 0.05')
legend('$\hat{a}_{1}$','$\hat{a}_{2}$','$a_{1}$','$a_{2}$','Interpreter','latex')
xlabel('Sample')
ylabel('Weight')
ylim([0 0.9])

subplot(L_u,1,2)
plot(w_ave(1,:,2),'LineWidth',1.5)
hold on
grid on
plot(w_ave(2,:,2),'LineWidth',1.5)
plot([0 n],[AR_coef(1) AR_coef(1)],'LineWidth',1.5,'LineStyle','--')
plot([0 n],[AR_coef(2) AR_coef(2)],'LineWidth',1.5,'LineStyle','--')
title('Updating Curve of the Adaptive Filter Coefficients for \mu = 0.01')
legend('$\hat{a}_{1}$','$\hat{a}_{2}$','$a_{1}$','$a_{2}$','Interpreter','latex')
xlabel('Sample')
ylabel('Weight')
ylim([0 0.9])

tilefigs([0 0.4 0.4 1])
%% Initialization
clc;
clear;
close all;

%% Setting Parameters
loop = 100;
n = 1000;  % number of samples
n_var = 0.5;  % noise variance
MA_coef = 0.9;  % MA coefficient
L_MA = length(MA_coef);
u0_gass = 0;  % initial learning rate
u0_gngd = 5;
rho = 0.01;  % step sizes

%% Benveniste GASS and GNGD
model = arima('Constant',0,'MA',MA_coef,'Variance',n_var);  % Create univariate autoregressive integrated moving average (ARIMA) model
[S_MA,E,~] = simulate(model,n,'NumPaths',loop);  % Simulate sample paths of the model
S_MA = S_MA';  % MA Signal
E = E';  % simulated innovation paths E
S_delay = zeros(L_MA+1,n);

% weight
w_gass = zeros(L_MA+1,n,loop);
w_gngd = zeros(L_MA+1,n,loop);
for i = 1:loop
    S = [0 S_MA(i,1:end-1)];
    for j = 1:L_MA+1
        S_delay(j,:) = [zeros(1,j) E(i,1:n-j)];
    end
    w_gass(:,:,i) = GASS(S_delay,S,u0_gass,rho,'Ben');
    w_gngd(:,:,i) = GNGD(S_delay,S,u0_gngd,rho);
end
w_gass_ave = mean(w_gass(:,:,:),3);
w_gngd_ave = mean(w_gngd(:,:,:),3);

%% Result
t = 0:n-1;

figure
plot(t,w_gass_ave(2,:),'LineWidth',1.5)
hold on
grid on
plot(t,w_gngd_ave(2,:),'LineWidth',1.5)
title('Weight Estimates')
legend('Benveniste’s GASS','GNGD')
xlabel('Sample')
ylabel('Magnitude')

figure
plot(t,w_gass_ave(2,:),'LineWidth',1.5)
hold on
grid on
plot(t,w_gngd_ave(2,:),'LineWidth',1.5)
title('Weight Estimates')
legend('Benveniste’s GASS','GNGD')
xlabel('Sample')
ylabel('Magnitude')
axis([0 50 0 1])

tilefigs([0 0.4 0.7 1])
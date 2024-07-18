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
u = [0.01 0.1];  % learning rate
u0 = 0;  % initial learning rate
rho = 0.001;  % step sizes

%% three GASS algorithms against Standard LMS
model = arima('Constant',0,'MA',MA_coef,'Variance',n_var);  % Create univariate autoregressive integrated moving average (ARIMA) model
[S_MA,E,~] = simulate(model,n,'NumPaths',loop);  % Simulate sample paths of the model

S_MA = S_MA';  % MA Signal
E = E';  % simulated innovation paths E
S_delay = zeros(L_MA+1,n);

% weight
w_LMS_1 = zeros(L_MA+1,n,loop);
w_LMS_2 = zeros(L_MA+1,n,loop);
w_Ben = zeros(L_MA+1,n,loop);
w_Ang = zeros(L_MA+1,n,loop);
w_Mat = zeros(L_MA+1,n,loop);

for i = 1:loop
    S = [0 S_MA(i,1:end-1)];
    for j = 1:L_MA+1
        S_delay(j,:) = [zeros(1,j) E(i,1:n-j)];
    end
    [w_LMS_1(:,:,i),~] = LMS_AR_1(S_delay,S,u(1));
    [w_LMS_2(:,:,i),~] = LMS_AR_1(S_delay,S,u(2));
    w_Ben(:,:,i) = GASS(S_delay,S,u0,rho,'Ben');
    w_Ang(:,:,i) = GASS(S_delay,S,u0,rho,'Ang');
    w_Mat(:,:,i) = GASS(S_delay,S,u0,rho,'Mat');
end
w_LMS_ave_1 = mean(w_LMS_1(:,:,:),3);
w_LMS_ave_2 = mean(w_LMS_2(:,:,:),3);
w_Ben_ave = mean(w_Ben(:,:,:),3);
w_Ang_ave = mean(w_Ang(:,:,:),3);
w_Mat_ave = mean(w_Mat(:,:,:),3);

%% Result
% u=0.01 u=0.1 B A M
transient_error = [pow2db(abs(mean((0.9-w_LMS_ave_1(2,800:end)).^2)));
                   pow2db(abs(mean((0.9-w_LMS_ave_2(2,800:end)).^2)));
                   pow2db(abs(mean((0.9-w_Ben_ave(2,800:end)).^2)));
                   pow2db(abs(mean((0.9-w_Ang_ave(2,800:end)).^2)));
                   pow2db(abs(mean((0.9-w_Mat_ave(2,800:end)).^2)))]

t = 0:n-1;

figure;
sgtitle('Weight Error Curves')

subplot(2,1,1)
plot(t,0.9-w_LMS_ave_1(2,:),'LineWidth',1.5)
hold on
grid on
plot(t,0.9-w_LMS_ave_2(2,:),'LineWidth',1.5)
plot(t,0.9-w_Ben_ave(2,:),'LineWidth',1.5)
plot(t,0.9-w_Ang_ave(2,:),'LineWidth',1.5)
plot(t,0.9-w_Mat_ave(2,:),'LineWidth',1.5)
title(['\rho = ' num2str(rho)])
legend('\mu = 0.01','\mu = 0.1','Benveniste','Ang & Farhang','Matthews & Xie')
ylabel('Magnitude')

subplot(2,1,2)
plot(t,0.9-w_LMS_ave_1(2,:),'LineWidth',1.5)
hold on
grid on
plot(t,0.9-w_LMS_ave_2(2,:),'LineWidth',1.5)
plot(t,0.9-w_Ben_ave(2,:),'LineWidth',1.5)
plot(t,0.9-w_Ang_ave(2,:),'LineWidth',1.5)
plot(t,0.9-w_Mat_ave(2,:),'LineWidth',1.5)
xlabel('Sample')
ylabel('Magnitude')
axis([0 100 0 1])

tilefigs([0 0.4 0.25 1])
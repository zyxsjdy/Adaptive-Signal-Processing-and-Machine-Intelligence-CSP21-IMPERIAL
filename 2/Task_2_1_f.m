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

% LMS leakage
leak = 0.1:0.1:0.9;
% number of leakages
L_L = length(leak);

%% LMS Adaptive Predictor
model = arima('Constant',0,'AR',AR_coef,'Variance',n_var);  % Create univariate autoregressive integrated moving average (ARIMA) model
S_AR = (simulate(model,n,'NumPaths',loop))';  % Simulate sample paths of the model

w = zeros(L_AR,n,loop);  % weight
w_ave = zeros(L_AR,n,L_L,L_u);  % average weight
for l = 1:L_L
    for k = 1:L_u
        for i = 1:loop
            S = S_AR(i,:);
            [w(:,:,i),~,~] = LMS_AR(S,L_AR,1,u(k),leak(l));
        end
        w_ave(:,:,l,k) = mean(w(:,:,:),3);
    end
end

%% Result
u5 = zeros(L_u,L_L);
u1 = zeros(L_u,L_L);
for l = 1:L_L
    u5(:,l) = [mean(w_ave(1,800:end,l,1)) mean(w_ave(2,800:end,l,1))]';
    u1(:,l) = [mean(w_ave(1,800:end,l,2)) mean(w_ave(2,800:end,l,2))]';
end
fprintf('coefficients with u = 0.05:\n')
fprintf([num2str(u5(1,:)) '\n'])
fprintf([num2str(u5(2,:)) '\n'])
fprintf('coefficients with u = 0.01:\n')
fprintf([num2str(u1(1,:)) '\n'])
fprintf([num2str(u1(2,:)) '\n'])


% figure
% sgtitle('Steady State Values of the Adaptive Filter Coefficients')
% 
% for l = 1:L_L
%     for k = 1:L_u
%         subplot(L_L,L_u,(l-1)*L_u+k)
%         plot(w_ave(1,:,l,k),'LineWidth',1.5)
%         hold on
%         grid on
%         plot(w_ave(2,:,l,k),'LineWidth',1.5)
%         plot([0 n],[AR_coef(1) AR_coef(1)],'LineWidth',1.5,'LineStyle','--')
%         plot([0 n],[AR_coef(2) AR_coef(2)],'LineWidth',1.5,'LineStyle','--')
%         title(['\mu = ' num2str(u(k)) ', \gamma = ' num2str(leak(l))])
%         legend('$\hat{a}_{1}$','$\hat{a}_{2}$','$a_{1}$','$a_{2}$','Interpreter','latex')
%         xlabel('Sample')
%         ylabel('Magnitude')
%         ylim([0 1])
%     end
% end
% 
% tilefigs([0 0.3 0.6 1])
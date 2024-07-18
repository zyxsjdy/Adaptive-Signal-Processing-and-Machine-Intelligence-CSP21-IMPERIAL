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
t = 0.5;
a = 1:t:200;  % scaling the activation function
n_a = length(a);

%% LMS using a*tanh
S_est = zeros(n_a,length(S));
MSPE = zeros(n_a,1);
R_p = zeros(n_a,1);
for i = 1:n_a
    [~,S_est(i,:),error] = LMS_tanh(S,p,delay,u,0,a(i));  % apply the LMS prediction
    MSPE(i) = pow2db(mean(abs(error).^2));  % MSPE
    R_p(i) = pow2db(var(S_est(i,:))/var(error));  % prediction gain
end
[M_MSPE,I_MSPE] = min(MSPE);
[M_R_p,I_R_p] = max(R_p);

%% Result
fprintf(['min MSPE = ' num2str(M_MSPE) ' dB when a = ' num2str((I_MSPE-1)*t+1) '\n'])
fprintf(['max R_p = ' num2str(M_R_p) ' dB when a = ' num2str((I_R_p-1)*t+1) '\n'])

a1 = [20 (I_MSPE-1)*t+1 100];
n_a1 = length(a1);
figure
sgtitle('Zero-mean Version of y and its One-step Ahead Prediction using LMS with a*tanh')
for i = 1:n_a1
    subplot(n_a1,1,i)
    plot(S,'Linewidth',1.5)
    hold on
    grid on
    plot(S_est(round((a1(i)-1)/t+1),:),'Linewidth',1.5)
    title(['a = ' num2str(a1(i))])
    legend('Zero-mean Version','One-step Ahead Prediction')
    ylabel('Magnitude')
end
xlabel('Sample')

figure
plot(a,MSPE,'LineWidth',1.5)
hold on
grid on
plot(a,R_p,'LineWidth',1.5)
title('MSPE and Prediction Gain of LMS with a*tanh')
legend('MSPE','R_p')
xlabel('Activation Function Scale')
ylabel('Magnitude (dB)')

tilefigs([0 0.4 0.7 1])
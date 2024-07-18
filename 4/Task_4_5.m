%% Initialization
clc;
clear;
close all;

%% Setting Parameters
load('time-series.mat');
S = y';
p = 4;  % AR order
delay = 1;
u = 1e-5;  % learning rate
t = 0.1;
a = 1:t:100;  % scaling the activation function
n_a = length(a);

%% biased LMS using a*tanh with pre-trained weights

n_S_s = 20;  % over-fitting to a small number of samples
epochs = 100;
S_s = S(1:n_S_s);

h0 = zeros(p+1,n_a);
for i = 1:n_a
    for j = 1:epochs
        [h,~,~] = LMS_tanh_bias(S_s,p,delay,u,0,a(i),h0(:,i));  % obtain the weight
        h0(:,i) = h(:,end);  % take the latest weight
    end
end


S_est = zeros(n_a,length(S));
MSPE = zeros(n_a,1);
MSPE_0 = zeros(n_a,1);
R_p = zeros(n_a,1);
R_p_0 = zeros(n_a,1);
for i = 1:n_a
    [~,S_est(i,:),error] = LMS_tanh_bias(S,p,delay,u,0,a(i),h0(:,i));  % apply the LMS prediction
    MSPE(i) = pow2db(mean(abs(error).^2));  % MSPE
    R_p(i) = pow2db(var(S_est(i,:))/var(error));  % prediction gain

    MSPE_0(i) = pow2db(mean(abs(error(1:200)).^2));  % MSPE in the first 200 samples
    R_p_0(i) = pow2db(var(S_est(i,1:200))/var(error(1:200)));  % prediction gain in the first 200 samples
end
[M_MSPE,I_MSPE] = min(MSPE);
[M_R_p,I_R_p] = max(R_p);

%% Result
fprintf(['min MSPE = ' num2str(M_MSPE) ' dB when a = ' num2str((I_MSPE-1)*t+1) '\n'])
fprintf(['max R_p = ' num2str(M_R_p) ' dB when a = ' num2str((I_R_p-1)*t+1) '\n'])
fprintf(['MSPE in the first 200 samples = ' num2str(MSPE_0(I_MSPE)) ' dB \n'])
fprintf(['R_p in the first 200 samples = ' num2str(R_p_0(I_MSPE)) ' dB \n'])

a1 = [20 (I_MSPE-1)*t+1 100];
n_a1 = length(a1);
figure
sgtitle('y and its One-step Ahead Prediction using LMS algorithm with activation function, bias and pre-trained weights')
for i = 1:n_a1
    subplot(n_a1,1,i)
    plot(S,'Linewidth',1.5)
    hold on
    grid on
    plot(S_est(round((a1(i)-1)/t+1),:),'Linewidth',1.5)
    title(['a = ' num2str(a1(i))])
    legend('y','One-step Ahead Prediction')
    ylabel('Magnitude')
end
xlabel('Sample')

figure
plot(a,MSPE,'LineWidth',1.5)
hold on
grid on
plot(a,R_p,'LineWidth',1.5)
title('MSPE and Prediction Gain of biased LMS using tanh and Pre-trained Weights')
legend('MSPE','R_p')
xlabel('Activation Function Scale')
ylabel('Magnitude (dB)')

tilefigs([0 0.4 0.7 1])
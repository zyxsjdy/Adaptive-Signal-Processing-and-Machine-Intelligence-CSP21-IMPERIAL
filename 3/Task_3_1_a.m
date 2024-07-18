%% Initialization
clc;
clear;
close all;

%% Setting Parameters
loop = 100;
n = 1000;  % number of samples
MA_coef = [1.5+1i 2.5-0.5i];  % MA coefficient
L_MA = length(MA_coef);
u = 0.1;  % step sizes

% signal generation
CWGN = 1/sqrt(2) * (randn(loop,n) + 1i*randn(loop,n));  % circular white Gaussian noise
WLMA = CWGN + MA_coef(1)*[zeros(loop,1) CWGN(:,1:end-1)] + MA_coef(2)*[zeros(loop,1) conj(CWGN(:,1:end-1))];  % WLMA(1)
fprintf(['circularity coefficient of WLMA: ' num2str(circularity(WLMA)) '\n'])

%% CLMS and ACLMS
CLMS_est = zeros(loop,n);
ACLMS_est = zeros(loop,n);
e_CLMS = zeros(loop,n);
e_ACLMS = zeros(loop,n);
CWGN_delay = zeros(L_MA,n);
for i = 1:loop
    for j = 1:L_MA
        CWGN_delay(j,:) = [zeros(1,j) CWGN(i,1:n-j)];  % delay the CWGN
    end
    WLMA_delay = [0 (WLMA(i,1:end-1))];  % delay the WLMA
    [CLMS_est(i,:),e_CLMS(i,:),~] = CLMS(CWGN_delay,WLMA_delay,u,0);  % apply CLMS
    [ACLMS_est(i,:),e_ACLMS(i,:),~,~] = ACLMS(CWGN_delay,WLMA_delay,u);  % apply ACLMS
end
e_CLMS_ave = pow2db(mean(abs(e_CLMS).^2,1));  % learning curve
e_ACLMS_ave = pow2db(mean(abs(e_ACLMS).^2,1));

%% Result
figure
subplot(1,2,1)
scatter(real(WLMA(:)),imag(WLMA(:)),2)
hold on
grid on
scatter(real(CLMS_est(:)),imag(CLMS_est(:)),0.2)
scatter(real(ACLMS_est(:)),imag(ACLMS_est(:)),0.2)
title('Identification of the WLMA(1) Model')
legend('WLMA(1)','CLMS','ACLMS')
xlabel('Real')
ylabel('Imag')
axis equal
axis([-15 15 -8 8])

subplot(1,2,2)
plot(e_CLMS_ave,'LineWidth',1.5)
hold on
grid on
plot(e_ACLMS_ave,'LineWidth',1.5)
title('Learning Curve for the CLMS and ACLMS')
legend('CLMS','ACLMS')
xlabel('Sample')
ylabel('Magnitude (dB)')
ylim([-350 50])

tilefigs([0 0.64 0.7 1])
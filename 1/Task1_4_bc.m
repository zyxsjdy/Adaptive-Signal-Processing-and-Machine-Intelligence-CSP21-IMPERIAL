%% Initialization
clc;
clear;
close all;

%% Setting Parameters
Fs = 1;  % normalized sampling frequency
n = 1000;  % number of samples

ARcoef = [2.76 -3.81 2.65 -0.92];  % coefficients
p = 2:14;  % model orders
l_p = length(p);

%% Realization (b)
x = zeros(1,n);
for i = 5:n
   x(i) = ARcoef(1)*x(i-1) + ARcoef(2)*x(i-2) + ARcoef(3)*x(i-3) + ARcoef(4)*x(i-4) + randn;    
end
x = x(500:end);
n = length(x);  % new number of samples

[h,f_b] = freqz(1,[1 -ARcoef],n,Fs);  % frequency response of digital filter
PSD_b = pow2db(abs(h).^2);  % True PSD

error_b = zeros(l_p,1);
PSD_est_b = zeros(l_p,n);
for i = 1:l_p
    [ARcoef_est,error_b(i)] = aryule(x,p(i));  % Yule-Walker Method
    [h,~] = freqz(sqrt(error_b(i)),ARcoef_est,n);
    PSD_est_b(i,:) = pow2db(abs(h).^2);
end

%% Realization (c)
n = 10000;  % number of samples
x = zeros(1,n);
for i = 5:n
   x(i) = ARcoef(1)*x(i-1) + ARcoef(2)*x(i-2) + ARcoef(3)*x(i-3) + ARcoef(4)*x(i-4) + randn;    
end
x = x(500:end);
n = length(x);  % new number of samples

[h,f_c] = freqz(1,[1 -ARcoef],n,Fs);  % frequency response of digital filter
PSD_c = pow2db(abs(h).^2);  % True PSD

error_c = zeros(l_p,1);
PSD_est_c = zeros(l_p,n);
for i = 1:l_p
    [ARcoef_est,error_c(i)] = aryule(x,p(i));  % Yule-Walker Method
    [h,~] = freqz(sqrt(error_c(i)),ARcoef_est,n);
    PSD_est_c(i,:) = pow2db(abs(h).^2);
end

%% Result
figure
sgtitle('PSD using AR model');

subplot(2,2,1)
plot(f_b,PSD_b,'LineWidth',1.5)
hold on
plot(f_b,PSD_est_b(1,:),'LineWidth',1.5)
grid on
legend('True PSD','Order 2')
xlabel('Normalised Frequency (× \pi rad/sample)')
ylabel('Magnitude (dB)')
subplot(2,2,2)
plot(f_b,PSD_b,'LineWidth',1.5)
hold on
plot(f_b,PSD_est_b(3,:),'LineWidth',1.5)
grid on
legend('True PSD','Order 4')
xlabel('Normalised Frequency (× \pi rad/sample)')
ylabel('Magnitude (dB)')
subplot(2,2,3)
plot(f_b,PSD_b,'LineWidth',1.5)
hold on
plot(f_b,PSD_est_b(5,:),'LineWidth',1.5)
grid on
legend('True PSD','Order 6')
xlabel('Normalised Frequency (× \pi rad/sample)')
ylabel('Magnitude (dB)')
subplot(2,2,4)
plot(f_b,PSD_b,'LineWidth',1.5)
hold on
plot(f_b,PSD_est_b(7,:),'LineWidth',1.5)
grid on
legend('True PSD','Order 8')
xlabel('Normalised Frequency (× \pi rad/sample)')
ylabel('Magnitude (dB)')

figure
plot(p, pow2db(error_b), 'LineWidth', 2)
grid on
legend('Error')
title('Error vs AR Order')
xlabel('Order')
ylabel('Magnitude (dB)')


figure
sgtitle('PSD using AR model');

subplot(2,2,1)
plot(f_c,PSD_c,'LineWidth',1.5)
hold on
plot(f_c,PSD_est_c(1,:),'LineWidth',1.5)
grid on
legend('True PSD','Order 2')
xlabel('Normalised Frequency (× \pi rad/sample)')
ylabel('Magnitude (dB)')
subplot(2,2,2)
plot(f_c,PSD_c,'LineWidth',1.5)
hold on
plot(f_c,PSD_est_c(3,:),'LineWidth',1.5)
grid on
legend('True PSD','Order 4')
xlabel('Normalised Frequency (× \pi rad/sample)')
ylabel('Magnitude (dB)')
subplot(2,2,3)
plot(f_c,PSD_c,'LineWidth',1.5)
hold on
plot(f_c,PSD_est_c(5,:),'LineWidth',1.5)
grid on
legend('True PSD','Order 6')
xlabel('Normalised Frequency (× \pi rad/sample)')
ylabel('Magnitude (dB)')
subplot(2,2,4)
plot(f_c,PSD_c,'LineWidth',1.5)
hold on
plot(f_c,PSD_est_c(7,:),'LineWidth',1.5)
grid on
legend('True PSD','Order 8')
xlabel('Normalised Frequency (× \pi rad/sample)')
ylabel('Magnitude (dB)')

figure
plot(p, pow2db(error_c), 'LineWidth', 2)
grid on
legend('Error')
title('Error vs AR Order')
xlabel('Order')
ylabel('Magnitude (dB)')


tilefigs([0 0.04 0.7 1])
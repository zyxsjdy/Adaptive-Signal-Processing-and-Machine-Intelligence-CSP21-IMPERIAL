%% Initialization
clc;
clear;
close all;

%% Setting Parameters
load('low-wind.mat');
v(1,:) = (v_east' + 1i*v_north');
load('medium-wind.mat');
v(2,:) = (v_east' + 1i*v_north');
load('high-wind.mat');
v(3,:) = (v_east' + 1i*v_north');
n_w = 3;
n = length(v(3,:));  % signal length
u = [1e-1 1e-2 1e-3];  % step sizes
M = 1:30;  % filter length
n_M = length(M);

%% CLMS and ACLMS
MSPE_CLMS = zeros(n_w,n_M);
MSPE_ACLMS = zeros(n_w,n_M);
for j = 1:n_w
    for i = 1:n_M
        S_delay = zeros(M(i),n);
        for k = 1:M(i)
            S_delay(k,:) = [zeros(1,k) v(j,1:n-k)];  % delay the wind data
        end

        [~,error,~] = CLMS(S_delay,v(j,:),u(j),0);  % apply CLMS
        MSPE_CLMS(j,i) = pow2db(mean(abs(error).^2));

        [~,error,~,~] = ACLMS(S_delay,v(j,:),u(j));  % apply ACLMS
        MSPE_ACLMS(j,i) = pow2db(mean(abs(error).^2));
    end
end

%% Result
figure
subplot(2,2,1)
scatter(real(v(1,:)),imag(v(1,:)),'.')
hold on
grid on
scatter(mean(real(v(1,:))),mean(imag(v(1,:))),100,'r','.') % show the data center
rho = circularity(v(1,:));
title(['Low-wind, |\rho| = ' num2str(round(rho,3))])
xlabel('Real')
ylabel('Imag')
axis equal
axis([-0.4 0.4 -0.4 0.4])
subplot(2,2,2)
scatter(real(v(2,:)),imag(v(2,:)),'.')
hold on
grid on
scatter(mean(real(v(2,:))),mean(imag(v(2,:))),100,'r','.') % show the data center
rho = circularity(v(2,:));
title(['Medium-wind, |\rho| = ' num2str(round(rho,3))])
xlabel('Real')
ylabel('Imag')
axis equal
axis([-2 2 -2 2])
subplot(2,2,3)
scatter(real(v(3,:)),imag(v(3,:)),'.')
hold on
grid on
scatter(mean(real(v(3,:))),mean(imag(v(3,:))),100,'r','.') % show the data center
rho = circularity(v(3,:));
title(['High-wind, |\rho| = ' num2str(round(rho,3))])
xlabel('Real')
ylabel('Imag')
axis equal
axis([-4 4 -4 4])
subplot(2,2,4)
scatter(real(v(3,:)),imag(v(3,:)),'b','.')
hold on
grid on
scatter(real(v(2,:)),imag(v(2,:)),'r','.')
scatter(real(v(1,:)),imag(v(1,:)),'y','.')
scatter(mean(real(v(3,:))),mean(imag(v(3,:))),100,'g','.')
scatter(mean(real(v(2,:))),mean(imag(v(2,:))),100,'g','.')
scatter(mean(real(v(1,:))),mean(imag(v(1,:))),100,'g','.')
title('Wind Data')
xlabel('Real')
ylabel('Imag')
axis equal
axis([-4 4 -4 4])

figure
sgtitle('MSPE versus Filter Length')
subplot(3,1,1)
plot(MSPE_CLMS(1,:),'LineWidth',1.5)
hold on
grid on
plot(MSPE_ACLMS(1,:),'LineWidth',1.5)
title('Low-wind')
legend('CLMS','ACLMS')
ylabel('Magnitude (dB)')
subplot(3,1,2)
plot(MSPE_CLMS(2,:),'LineWidth',1.5)
hold on
grid on
plot(MSPE_ACLMS(2,:),'LineWidth',1.5)
title('Medium-wind')
legend('CLMS','ACLMS')
ylabel('Magnitude (dB)')
subplot(3,1,3)
plot(MSPE_CLMS(3,:),'LineWidth',1.5)
hold on
grid on
plot(MSPE_ACLMS(3,:),'LineWidth',1.5)
title('High-wind')
legend('CLMS','ACLMS')
xlabel('Filter Length')
ylabel('Magnitude (dB)')


tilefigs([0 0.4 0.7 1])
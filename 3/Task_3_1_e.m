%% Initialization
clc;
clear;
close all;

%% Setting Parameters
Fs = 1000;  % sampling frequency
n = 1000;  % number of samples
t = (0:n-1)/Fs;  % sampling time

f0 = 50;  % system frequency
AMP0 = [1.0 1.0 1.0];
phi = 0;  % phase shift
delta0 =[0 0];  % phase distortions
u = 0.05;  % step size

%% Balanced System
Va = AMP0(1) * cos(2*pi*f0*t + phi);
Vb = AMP0(2) * cos(2*pi*f0*t + phi + delta0(1) - (2*pi/3));
Vc = AMP0(3) * cos(2*pi*f0*t + phi + delta0(2) + (2*pi/3));
V3_b = [Va;Vb;Vc];  % three-phase voltages vector
V3_b_CT = CLARKE(V3_b);  % Clarke Transform
V_b = V3_b_CT(2,:) + 1i*V3_b_CT(3,:);  % complex Clarke voltage under balanced conditions

V_b_delay = [0 V_b(1:n-1)];  % delay the signal
[~,~,h] = CLMS_1(V_b_delay,V_b,u,0);  % apply CLMS
f_CLMS_b = abs(Fs/(2*pi)*atan(imag(h)./real(h)));
[~,~,h,g] = ACLMS_1(V_b_delay,V_b,u);  % apply ACLMS
f_ACLMS_b = abs(Fs/(2*pi)*atan(sqrt(imag(h).^2-abs(g).^2)./real(h)));

%% Unbalanced System
% unbalanced amplitude
AMP = [0.2  1  1.8];  
V3_ubA = AMP'.*V3_b;  % three-phase voltages vector, unbalanced amplitude
V3_ubA_CT = CLARKE(V3_ubA);  % Clarke Transform
V_ubA = V3_ubA_CT(2,:) + 1i*V3_ubA_CT(3,:);  % complex Clarke voltage
Circu_V_ubA = circularity(V_ubA);  % circularity

V_ubA_delay = [0 V_ubA(1:n-1)];  % delay the signal
[~,~,h] = CLMS_1(V_ubA_delay,V_ubA,u,0);  % apply CLMS
f_CLMS_ubA = abs(Fs/(2*pi)*atan(imag(h)./real(h)));
[~,~,h,g] = ACLMS_1(V_ubA_delay,V_ubA,u);  % apply ACLMS
f_ACLMS_ubA = abs(Fs/(2*pi)*atan(sqrt(imag(h).^2-abs(g).^2)./real(h)));


% unbalanced phase
delta = [-pi/4  pi/4];  
Va_ubP = AMP0(1) * cos(2*pi*f0*t + phi);
Vb_ubP = AMP0(2) * cos(2*pi*f0*t + phi + delta(1) - (2*pi/3));
Vc_ubP = AMP0(3) * cos(2*pi*f0*t + phi + delta(2) + (2*pi/3));
V3_ubP = [Va_ubP;Vb_ubP;Vc_ubP];  % three-phase voltages vector, unbalanced phase
V3_ubP_CT = CLARKE(V3_ubP);  % Clarke Transform
V_ubP = V3_ubP_CT(2,:) + 1i*V3_ubP_CT(3,:);  % complex Clarke voltage
Circu_V_ubP = circularity(V_ubP);  % circularity

V_ubP_delay = [0 V_ubP(1:n-1)];  % delay the signal
[~,~,h] = CLMS_1(V_ubP_delay,V_ubP,u,0);  % apply CLMS
f_CLMS_ubP = abs(Fs/(2*pi)*atan(imag(h)./real(h)));
[~,~,h,g] = ACLMS_1(V_ubP_delay,V_ubP,u);  % apply ACLMS
f_ACLMS_ubP = abs(Fs/(2*pi)*atan(sqrt(imag(h).^2-abs(g).^2)./real(h)));

%% Result
figure
subplot(3,1,1)
plot(f_CLMS_b,'LineWidth',1.5)
hold on
grid on
plot(f_ACLMS_b,'LineWidth',1.5)
yline(50,'LineWidth',1.5,'Color','k','LineStyle','--')
title('Frequency Estimate for the Balanced System, |\rho| = 0.000')
legend('CLMS','ACLMS','f_0')
xlabel('Sample')
ylabel('frequency (Hz)')
ylim([0 60])

subplot(3,1,2)
plot(f_CLMS_ubA,'LineWidth',1.5)
hold on
grid on
plot(f_ACLMS_ubA,'LineWidth',1.5)
yline(50,'LineWidth',1.5,'Color','k','LineStyle','--')
title(['Frequency Estimate for the System with Unbalanced Amplitude, |\rho| = ' num2str(Circu_V_ubA,3)])
legend('CLMS','ACLMS','f_0')
xlabel('Sample')
ylabel('frequency (Hz)')
ylim([0 60])

subplot(3, 1, 3)
plot(f_CLMS_ubP,'LineWidth',1.5)
hold on
grid on
plot(f_ACLMS_ubP,'LineWidth',1.5)
yline(50,'LineWidth',1.5,'Color','k','LineStyle','--')
title(['Frequency Estimate for the System with Unbalanced Phase, |\rho| = ' num2str(Circu_V_ubP,3)])
legend('CLMS','ACLMS','f_0')
xlabel('Sample')
ylabel('frequency (Hz)')
ylim([0 60])

tilefigs([0 0.4 0.4 1])
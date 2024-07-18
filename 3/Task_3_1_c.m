%% Initialization
clc;
clear;
close all;

%% Setting Parameters
Fs = 10000;  % sampling frequency
n = 1000;  % number of samples
t = (0:n-1)/Fs;  % sampling time

f0 = 50;  % system frequency
AMP0 = [1 1 1];
phi = 0;  % phase shift
delta0 =[0 0];  % phase distortions

%% Balanced System
Va = AMP0(1) * cos(2*pi*f0*t + phi);
Vb = AMP0(2) * cos(2*pi*f0*t + phi + delta0(1) - (2*pi/3));
Vc = AMP0(3) * cos(2*pi*f0*t + phi + delta0(2) + (2*pi/3));
V3_b = [Va;Vb;Vc];  % three-phase voltages vector
V3_b_CT = CLARKE(V3_b);  % Clarke Transform
V_b = V3_b_CT(2,:) + 1i*V3_b_CT(3,:);  % complex Clarke voltage under balanced conditions

%% Unbalanced System
AMP = [0.5  1    1.5;  % unbalanced amplitude
       0.5  1.5  1; 
       1.5  1    0.5; 
       0.2  1    1.8];

delta = [-pi/6  pi/6;  % unbalanced phase
          pi/6 -pi/6; 
          pi/6  pi/3; 
         -pi/4  pi/4];

L = 4;  % 4 trials
V_ubA = zeros(L,n);  % complex Clarke voltage with unbalanced amplitude
V_ubP = zeros(L,n);  % complex Clarke voltage with unbalanced phase

Circu_V_ubA = zeros(L,1);  % circularity
Circu_V_ubP = zeros(L,1);

for i = 1:L
    V3_ubA = AMP(i,:)'.*V3_b;  % three-phase voltages vector, unbalanced amplitude
    V3_ubA_CT = CLARKE(V3_ubA);  % Clarke Transform
    V_ubA(i,:) = V3_ubA_CT(2,:) + 1i*V3_ubA_CT(3,:);  % complex Clarke voltage
    Circu_V_ubA(i) = circularity(V_ubA(i,:));  % circularity

    Va_ubP = AMP0(1) * cos(2*pi*f0*t + phi);
    Vb_ubP = AMP0(2) * cos(2*pi*f0*t + phi + delta(i,1) - (2*pi/3));
    Vc_ubP = AMP0(3) * cos(2*pi*f0*t + phi + delta(i,2) + (2*pi/3));
    V3_ubP = [Va_ubP;Vb_ubP;Vc_ubP];  % three-phase voltages vector, unbalanced phase
    V3_ubP_CT = CLARKE(V3_ubP);  % Clarke Transform
    V_ubP(i,:) = V3_ubP_CT(2,:) + 1i*V3_ubP_CT(3,:);  % complex Clarke voltage
    Circu_V_ubP(i) = circularity(V_ubP(i,:));  % circularity
end

%% Result
figure

subplot(1,2,1)
scatter(real(V_b),imag(V_b),100,'.')
hold on
grid on
for i = 1:L
    scatter(real(V_ubA(i,:)),imag(V_ubA(i,:)),100,'.')
end
title('Circularity Diagram with Unbalanced Magnitude')
legend(['V_a = 1, V_b = 1, V_c = 1, |\rho| = 0.000    '; ...
        'V_a = 0.5, V_b = 1, V_c = 1.5, |\rho| = ' num2str(Circu_V_ubA(1),3); ...
        'V_a = 0.5, V_b = 1.5, V_c = 1, |\rho| = ' num2str(Circu_V_ubA(2),3); ...
        'V_a = 1.5, V_b = 1, V_c = 0.5, |\rho| = ' num2str(Circu_V_ubA(3),3); ...
        'V_a = 0.2, V_b = 1, V_c = 1.8, |\rho| = ' num2str(Circu_V_ubA(4),3)])
xlabel('Real')
ylabel('Imag')
axis equal
axis([-2 2 -2 2])


subplot(1,2,2)
scatter(real(V_b),imag(V_b),100,'.')
hold on
grid on
for i = 1:L
    scatter(real(V_ubP(i,:)),imag(V_ubP(i,:)),100,'.')
end
title('Circularity Diagram with Unbalanced Phases')
legend(['\Delta_b = 0, \Delta_c = 0, |\rho| = 0.000         '; ...
        '\Delta_b = -\pi/6, \Delta_c = \pi/6, |\rho| = ' num2str(Circu_V_ubP(1),3); ...
        '\Delta_b = \pi/6, \Delta_c = -\pi/6, |\rho| = ' num2str(Circu_V_ubP(2),3) '  '; ...
        '\Delta_b = \pi/6, \Delta_c = \pi/3, |\rho| = '  num2str(Circu_V_ubP(3),3) '   '; ...
        '\Delta_b = -\pi/4, \Delta_c = \pi/4, |\rho| = ' num2str(Circu_V_ubP(4),3)])
xlabel('Real')
ylabel('Imag')
axis equal
axis([-2 2 -2 2])

tilefigs([0 0.3 0.8 1])
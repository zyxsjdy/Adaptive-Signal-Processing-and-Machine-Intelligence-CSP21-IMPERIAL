%% Initialization
clc;
clear;
close all;

%% Setting Parameters
load('EEG_Data_Assignment1.mat');
n = 1200;  % number of samples
start = 1000;
POz = POz(start:start+n-1);
u = 1;  % step size

%% DFT-CLMS
% complex phasor
N_DFT = 1024;
x = zeros(N_DFT,n);
for k = 1:n
    for j = 1:N_DFT
        x(j,k) = 1/N_DFT * exp(1i*2*k*(j-1)*pi/N_DFT);
    end
end

[~,~,h] = CLMS(x,POz,u,0);  % apply CLMS
H = abs(h).^2;
medianH = 500*median(median(H));  % Remove outliers in the matrix H
H(H>medianH) = medianH;

%% Result
f1 = (0:N_DFT-1)/N_DFT * fs;

figure;
surf((1:n)/fs,f1,H,'LineStyle','none')
view(2);
cb = colorbar;
cb.Label.String = 'PSD';
title('Time Frequency Spectrum of the EEG Signal')
xlabel('Time (s)');
ylabel('frequency (Hz)');
ylim([0 60]);

tilefigs([0 0.5 0.6 1])
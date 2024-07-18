%% Initialization
clc;
clear;
close all;

%% Setting Parameters
n = 1500;  % number of samples
Fs = 1500;  % sampling frequency
M = 1;  % filter length
u = 0.1;  % step size

%% Use CLMS Algorithm to Estimate the AR Coefficient
% Generate the frequency function
f = 100*ones(1,1500);
for i = 501:1500
    if i <= 1000
        f(i) = f(i) + (i-500)/2;
    else
        f(i) = f(i) + ((i-1000)/25)^2;
    end
end
phi = cumsum(f);  % find the phase
y = exp(1i*2*pi/Fs*phi) + sqrt(0.05/2)*(randn(1,n)+1i*randn(1,n));  % FM signal
fprintf(['circularity coefficient of y: ' num2str(circularity(y)) '\n'])

y_delay = [0 y(1:n-1)];  % delay the signal
[~,~,h] = CLMS(y_delay,y,u,0);  % apply CLMS to estimate coefficients

N_DFT = 1024;
H = zeros(N_DFT,n);
for i = 1:n
    [h1,f1] = freqz(1,[1;-conj(h(i))],N_DFT,Fs);
    H(:,i) = abs(h1).^2;
end
medianH = 50*median(median(H));  % Remove outliers in the matrix H
H(H>medianH) = medianH;

%% Result
% Plot time-frequency diagram
surf(1:n,f1,H,'LineStyle','none')
view(2)
cb = colorbar;
cb.Label.String = 'PSD';
title(['Time Frequency Diagram using CLMS Algorithm with \mu = ' num2str(u)])
xlabel('Sample')
ylabel('frequency (Hz)')
ylim([0 600])
hold on
xf = 1:n;
yf = f;
zf = 100*ones(1,n);
plot3(xf,yf,zf,'Color','r','LineWidth',1.5)

tilefigs([0 0.5 0.6 1])
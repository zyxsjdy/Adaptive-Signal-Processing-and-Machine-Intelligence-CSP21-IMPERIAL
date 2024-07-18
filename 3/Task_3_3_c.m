%% Initialization
clc;
clear;
close all;

%% Setting Parameters
n = 1500;  % number of samples
Fs = 1500;  % sampling frequency
u = 1;  % step size
leak0 = 0;  % LMS leakage
leak1 = 0.1;

%% DFT-CLMS
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

% complex phasor
N_DFT = 1024;
x = zeros(N_DFT,n);
for k = 1:n
    for j = 1:N_DFT
        x(j,k) = 1/N_DFT * exp(1i*2*k*(j-1)*pi/N_DFT);
    end
end

[~,~,h] = CLMS(x,y,u,leak0);  % apply CLMS
H0 = abs(h).^2;
medianH = 50*median(median(H0));  % Remove outliers in the matrix H
H0(H0>medianH) = medianH;

[~,~,h] = CLMS(x,y,u,leak1);  % apply CLMS with leak
H1 = abs(h).^2;
medianH = 50*median(median(H1));  % Remove outliers in the matrix H
H1(H1>medianH) = medianH;

%% Result plot
f1 = (0:N_DFT-1)/N_DFT * Fs;

% Plot time-frequency diagram
figure
surf(1:n,f1,H0,'LineStyle','none')
view(2)
cb = colorbar;
cb.Label.String = 'PSD';
title(['Time Frequency Diagram using DFT-CLMS Algorithm with \mu = ' num2str(u)])
xlabel('Sample')
ylabel('frequency (Hz)')
ylim([0 600])
hold on
xf = 1:n;
yf = f;
zf = 100*ones(1,n);
plot3(xf,yf,zf,'Color','r','LineWidth',1.5)

figure
surf(1:n,f1,H1,'LineStyle','none')
view(2)
cb = colorbar;
cb.Label.String = 'PSD';
title(['Time Frequency Diagram using DFT-CLMS Algorithm with \mu = ' num2str(u) ' and \gamma = ' num2str(leak1)])
xlabel('Sample')
ylabel('frequency (Hz)')
ylim([0 600])
hold on
xf = 1:n;
yf = f;
zf = 100*ones(1,n);
plot3(xf,yf,zf,'Color','r','LineWidth',1.5)


tilefigs([0 0.5 1 1])
%% Initialization
clc;
clear;
close all;

%% Low-rank Approximation
load('PCAPCR.mat');

[U,S,V] = svd(Xnoise);  % singular value of Xnoise
rank_X = rank(X);  % rank of X

Xnoise_bar = U(:,1:rank_X) * S(1:rank_X,1:rank_X) * V(:,1:rank_X)';  % low-rank approximation

error_Xn = abs(vecnorm(X-Xnoise)).^2;
error_Xnbar = abs(vecnorm(X-Xnoise_bar)).^2;

%% Result
figure
stem(error_Xn)
hold on
grid on
stem(error_Xnbar)
title('Difference between the Variables')
legend('Noise Corrupted Matrix','Denoised Matrix')
legend('$X$ and $X_{noise}$','$X$ and $\widetilde{X}_{noise}$','Interpreter','latex')
xlabel('Index')
ylabel('Square Error')

tilefigs([0 0.4 0.4 1])
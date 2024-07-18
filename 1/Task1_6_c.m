%% Initialization
clc;
clear;
close all;
load('PCAPCR.mat');

%% Ordinary Least Squares (OLS)
B_OLS = (Xnoise'*Xnoise) \ Xnoise'*Y;  % OLS estimate for the unknown regression matrix

% output data
Y_OLS = Xnoise*B_OLS;
error_Y_OLS = abs(vecnorm(Y-Y_OLS)).^2;

Ytest_OLS = Xtest*B_OLS;
error_Ytest_OLS = abs(vecnorm(Ytest-Ytest_OLS)).^2;

%% Principle Component Regression (PCR)
rank_X = rank(X);  % rank of X

[U_n,S_n,V_n] = svd(Xnoise);  % singular value of Xnoise
Xnoise_bar = U_n(:,1:rank_X) * S_n(1:rank_X,1:rank_X) * V_n(:,1:rank_X)';  % low-rank approximation

B_PCR = V_n(:,1:rank_X) / S_n(1:rank_X,1:rank_X) * U_n(:,1:rank_X)' * Y;  % PCR solution

Y_PCR = Xnoise_bar*B_PCR;  % estimated data
error_Y_PCR = abs(vecnorm(Y-Y_PCR)).^2;  % error

[U_t,S_t,V_t] = svd(Xtest);  % singular value of Xtest
Xtest_bar = U_t(:,1:rank_X) * S_t(1:rank_X,1:rank_X) * V_t(:,1:rank_X)';

Ytest_PCR = Xtest_bar*B_PCR;  % estimated data
error_Ytest_PCR = abs(vecnorm(Ytest-Ytest_PCR)).^2;  % error

%% Result
fprintf(['error OLS using Y: ' num2str(error_Y_OLS,6) '\n'])
fprintf(['error PCR using Y: ' num2str(error_Y_PCR,6) '\n'])
fprintf(['error OLS using Y_test: ' num2str(error_Ytest_OLS,6) '\n'])
fprintf(['error PCR using Y_test: ' num2str(error_Ytest_PCR,6) '\n'])

%% Initialization
clc;
clear;
close all;

%% Setting Parameters
load('PCAPCR.mat');
loop = 1000;

B_OLS = (Xnoise'*Xnoise) \ Xnoise'*Y;  % OLS estimate for the unknown regression matrix

rank_X = rank(X);  % rank of X
[U,S,V] = svd(Xnoise);  % singular value of Xnoise
B_PCR = V(:,1:rank_X) / S(1:rank_X,1:rank_X) * U(:,1:rank_X)' * Y;  % PCR solution

%% Compute MSE
error_OLS = zeros(loop,5);
error_PCR = zeros(loop,5);
for i = 1:loop
    [Ytest,Ytest_OLS] = regval(B_OLS);  % estimated data
    error_OLS(i,:) = abs(vecnorm(Ytest-Ytest_OLS)).^2;  % error

    [Ytest,Ytest_PCR] = regval(B_PCR);  % estimated data
    error_PCR(i,:) = abs(vecnorm(Ytest-Ytest_PCR)).^2;  % error
end
errorOlsAvg = mean(error_OLS);
errorPcrAvg = mean(error_PCR);

%% Result
fprintf(['error OLS: ' num2str(errorOlsAvg,6) '\n'])
fprintf(['error PCR: ' num2str(errorPcrAvg,6) '\n'])

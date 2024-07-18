%% Initialization
clc;
clear;
close all;

%% SVD
load('PCAPCR.mat');

SVD_X = svd(X);  % singular value of X
SVD_X_n = svd(Xnoise);  % singular value of Xnoise

rank_X = rank(X);  % rank of X
rank_X_n = rank(Xnoise);  % rank of Xnoise

%% Result
fprintf(['Rank of X is ' num2str(rank_X) '\n'])
fprintf(['Rank of X_noise is ' num2str(rank_X_n) '\n'])


figure
stem(SVD_X_n)
hold on
grid on
stem(SVD_X)
legend('Xnoise','X')
title('Singular Value of X and Xnoise')
xlabel('Index of Singular Value')
ylabel('Magnitude')

tilefigs([0 0.4 0.4 1])
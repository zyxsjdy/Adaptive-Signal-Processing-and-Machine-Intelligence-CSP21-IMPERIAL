function [w,S_est] = LMS_AR_1(S_delay,S,u)

[p,n] = size(S_delay);  % take the size of the input
w = zeros(p,n+1);  % initialize the weight
S_est = zeros(1,n);  % initialize the estimated signal
error = zeros(1,n);  % initialize the error

for i = 1:n
    S_est(i) = w(:,i)' * S_delay(:,i);  % multiply the weight and input to get the estimated signal
    error(i) = S(i) - S_est(i);  % find the error between estimated signal and the desired signal
    w(:,i+1) = w(:,i) + S_delay(:,i)*error(i)*u;  % update the weight
end
w = w(:,2:end);  % get rid of the first 0
end
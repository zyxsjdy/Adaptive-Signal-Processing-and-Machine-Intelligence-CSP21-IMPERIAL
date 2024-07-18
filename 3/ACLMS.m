function [S_est,error,h,g] = ACLMS(S_delay,S,u)
[L,n] = size(S_delay);  % take the size of the input
h = zeros(L,n+1);  % initialize the weight
g = zeros(L,n+1);  % initialize the weight
S_est = zeros(1,n);  % initialize the estimated signal
error = zeros(1,n);  % initialize the error

for i = 1:n
    S_est(i) = h(:,i)'*S_delay(:,i) + g(:,i)'*conj(S_delay(:,i));  % multiply the weight and input to get the estimated signal
    error(i) = S(i) - S_est(i);  % find the error between estimated signal and the desired signal
    h(:,i+1) = h(:,i) + u*conj(error(i))*S_delay(:,i);  % update the weight
    g(:,i+1) = g(:,i) + u*conj(error(i))*conj(S_delay(:,i));  % update the weight
end
h = h(:,2:end);  % get rid of the first 0
g = g(:,2:end);
end
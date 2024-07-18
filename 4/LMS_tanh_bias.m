function [w,S_est,error] = LMS_tanh_bias(S,p,delay,u,leak,a,w)

n = length(S);
S_delay = zeros(p,n);
for i = 1:p
    S_delay(i,:) = [zeros(1,i+delay-1) S(1:n-(i+delay-1))];  % delay the signal
end
S_delay_bias = [ones(1,n);S_delay];  % augmented input

w(p+1,n+1) = 0;  % extend the weight
S_est = zeros(1,n);  % initialize the estimated signal
error = zeros(1,n);  % initialize the error

for i = 1:n
    S_est(i) = a*tanh(w(:,i)'*S_delay_bias(:,i));  % multiply the weight and input to get the estimated signal
    error(i) = S(i) - S_est(i);  % find the error between estimated signal and the desired signal
    w(:,i+1) = (1-leak*u)*w(:,i) + S_delay_bias(:,i)*error(i)*u;  % update the weight
end
w = w(:,2:end);  % get rid of the first 0
end
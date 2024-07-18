function w = GASS(S_delay,S,u,rho,alg)

[p,n] = size(S_delay);  % take the size of the input
w = zeros(p,n+1);  % initialize the weight
S_est = zeros(1,n);  % initialize the estimated signal
Error = zeros(1,n);  % initialize the error
J = zeros(p,n+1);  % initialize the cost function

for i = 1:n
    S_est(i) = w(:,i)' * S_delay(:,i);  % multiply the weight and input to get the estimated signal
    Error(i) = S(i) - S_est(i);  % find the error between estimated signal and the desired signal
    w(:,i+1) = w(:,i) + S_delay(:,i)*Error(i)*u;  % update the weight
    
    % Calculate the cost function
    if alg == 'Ben'
        J(:,i+1) = (eye(p) - S_delay(:,i)*S_delay(:,i)'*u) * J(:,i) + S_delay(:,i)*Error(i);
    elseif alg == 'Ang'
        J(:,i+1) = 0.6*J(:,i) + S_delay(:,i)*Error(i);
    elseif alg == 'Mat'
        J(:,i+1) = S_delay(:,i)*Error(i);
    else
        error('error in GASS');
    end
    u = u + rho*Error(i)*S_delay(:,i)'*J(:,i+1);  % update the step size
end
w = w(:,2:end);  % get rid of the first 0
end
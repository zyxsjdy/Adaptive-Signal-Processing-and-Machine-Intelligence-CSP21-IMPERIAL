function rho = circularity(signal)
COV = mean(abs(signal).^2);  % calculate the covariance
P_COV = mean(signal.^2);  % calculate the pseudocovariance

rho = abs(P_COV) / COV;  % circularity coefficient
end
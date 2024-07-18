function [V3_CT] = CLARKE(V3)
Clarke_Matrix = sqrt(2/3) * [1/sqrt(2)  1/sqrt(2)  1/sqrt(2);   % generate the CLARKE matrix
                             1          -1/2       -1/2; 
                             0          sqrt(3)/2  -sqrt(3)/2];
V3_CT = Clarke_Matrix * V3;  % apply the CLARKE transform
end
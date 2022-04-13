function [H] = exponential_10_Hess(x)
%EXPONENTIAL_10_HESS Summary of this function goes here
%   Detailed explanation goes here
H = zeros(10,10);
H(1,1) = 2*(exp(x(1)) * (exp(x(1)) + 1)^2 - (exp(x(1)) + 1) * exp(x(1)))  / (exp(x(1)) + 1)^4 + 0.1*exp(-x(1));
       %H(1,1) = 2*(-exp(2*x(1))+exp(x(1)))/(exp(x(1))+1)^3 + 0.1*exp(-x(1));
for i = 2:10
    H(i,i) = 12*(x(i) - 1)^2;
end
end


function [H] = exponential_1000_Hess(x)
%EXPONENTIAL_1000_HESS Summary of this function goes here
%   Detailed explanation goes here

H = diag(12*(x-1).^2);
H(1,1) = 2*(-exp(2*x(1))+exp(x(1)))/(exp(x(1))+1)^3 + 0.1*exp(-x(1));

end


function [g] = exponential_1000_grad(x)
%EXPONENTIAL_1000_GRAD Summary of this function goes here
%   Detailed explanation goes here

g = 4*(x-1).^3;
g(1) = 2*exp(x(1)) / (exp(x(1)) + 1)^2 - 0.1*exp(-x(1));

end


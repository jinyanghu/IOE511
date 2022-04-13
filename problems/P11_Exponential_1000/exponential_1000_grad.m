function [g] = exponential_1000_grad(x)
%EXPONENTIAL_1000_GRAD Summary of this function goes here
%   Detailed explanation goes here
g = zeros(1000,1);
g(1) = 2*exp(x(1)) / (exp(x(1)) + 1)^2 - 0.1*exp(-x(1));
for i = 2:1000
    g(i) = 4*(x(i)-1)^3;
end
end


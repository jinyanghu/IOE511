function [g] = exponential_10_grad(x)
%EXPONENTIAL_10_GRAD Summary of this function goes here
%   Detailed explanation goes here
g = zeros(10,1);
g(1) = 2*exp(x(1)) / (exp(x(1)) + 1)^2 - 0.1*exp(-x(1));
for i = 2:10
    g(i) = 4*(x(i)-1)^3;
end
end


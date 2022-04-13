function [f] = exponential_1000_func(x)
%EXPONENTIAL_1000_FUNC Summary of this function goes here
%   Detailed explanation goes here
z1 = (exp(x(1)) - 1)/(exp(x(1)) + 1) + 0.1 * exp(-x(1));
rest_sum = 0;
for i = 2:1000
    rest_sum = rest_sum + (x(i) - 1)^4;
end
f = z1 + rest_sum;
end


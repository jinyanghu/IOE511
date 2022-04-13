function [f] = rosenbrock_100_func(x)
%ROSENBROCK_100_FUNC Summary of this function goes here
%   Detailed explanation goes here
f = 0;
for i = 1:99
    f = f + (1 - x(i))^2 + 100*(x(i+1) - x(i)^2)^2;
end
end


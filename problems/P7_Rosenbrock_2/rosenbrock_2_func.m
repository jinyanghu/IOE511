function [f] = rosenbrock_2_func(x)
%ROSENBROCK_2_FUNC Summary of this function goes here
%   Detailed explanation goes here
f = (1 - x(1))^2 + 100*(x(2) - x(1)^2)^2;
end


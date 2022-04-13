function [g] = rosenbrock_100_grad(x)
%ROSENBROCK_100_GRAD Summary of this function goes here
%   Detailed explanation goes here
g = zeros(100,1);
g(1) = -2*(1-x(1)) - 400*(x(2)-x(1)^2)*x(1);
for i = 2:99
    g(i) = -2*(1-x(i)) - 400*(x(i+1)-x(i)^2)*x(i);
    g(i) = g(i) + 200 * (x(i) - x(i-1)^2);
end
g(100) = 200 * (x(100) - x(99)^2);
end


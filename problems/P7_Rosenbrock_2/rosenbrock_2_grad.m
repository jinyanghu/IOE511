function [g] = rosenbrock_2_grad(x)
%ROSENBROCK_2_GRAD Summary of this function goes here
%   Detailed explanation goes here
g_1 = -2*(1-x(1)) - 400*(x(2)-x(1)^2)*x(1);
g_2 = 200 * (x(2) - x(1)^2);
g = [g_1;g_2];
end


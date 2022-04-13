function [H] = rosenbrock_2_Hess(x)
%ROSENBROCK_2_HESS Summary of this function goes here
%   Detailed explanation goes here
h11 = 2-400*x(2) + 1200 * x(1)^2;
h22 = 200;
h12 = -400*x(1);
H = [h11 h12;h12 h22];
end


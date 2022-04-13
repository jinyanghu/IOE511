function [H] = rosenbrock_100_Hess(x)
%ROSENBROCK_100_HESS Summary of this function goes here
%   Detailed explanation goes here
H = zeros(100,100);
H(1,1) = 2-400*x(2) + 1200 * x(1)^2;
H(1,2) = -400*x(1);
H(2,1) = -400*x(1);
for i = 2:99
    %h11 = h11 + 2-400*x(2) + 1200 * x(1)^2;
    %h22 = h22+ 200;
    %h12 = -400*x(1);
    H(i,i) = 2-400*x(i) + 1200 * x(i)^2 + 200;
    H(i,i+1) = -400*x(i);
    H(i+1,i) = -400*x(i);
end
H(100,100) = 200;
end


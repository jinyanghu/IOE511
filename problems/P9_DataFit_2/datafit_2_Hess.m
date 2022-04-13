function [H] = datafit_2_Hess(x)
%DATAFIT_2_HESS Summary of this function goes here
%   Detailed explanation goes here
y = [1.5;2.25;2.625];
h11 = 0;
h12 = 0;
h22 = 0;
for i = 1:3
    h11 = h11+ 2*(y(i) - x(1)*(1-x(2)^i))*(1-x(2)^i)^2;
    h12 = h12 + 2*i*x(2)^(i-1)*(y(i)- x(1)*(1-x(2)^i)-x(1)*((1-x(2)^i)));
    h22 = h22 + 2*x(1)*i*((y(i) - x(1)*(1-x(2)^i))*(i-1)*x(2)^(i-2)+x(1)*i*x(2)^(2*i-2));
end    
H = [h11 h12;h12 h22];
end


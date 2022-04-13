function [f] = datafit_2_func(x)
%DATAFIT_2_FUNC Summary of this function goes here
y = [1.5;2.25;2.625];
f = 0;
for i = 1:3
    f = f + (y(i) - x(1)*(1-x(2)^i))^2;
end
end


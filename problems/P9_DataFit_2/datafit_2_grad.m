function [g] = datafit_2_grad(x)
y = [1.5;2.25;2.625];
g1 = 0;
g2 = 0;
for i= 1:3
   g1 = g1 + 2*(y(i) - x(1)*(1- x(2)^i))*(-(1-x(2)^i));
   g2 = g2 + 2*(y(i) - x(1)*(1- x(2)^i))*x(1)*i*x(2)^(i-1);
end
g = [g1;g2];
end


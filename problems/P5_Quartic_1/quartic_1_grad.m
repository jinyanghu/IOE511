function [g] = quartic_1_grad(x)
% Matrix Q
Q = [5 1 0 0.5;
     1 4 0.5 0;
     0 0.5 3 0;
     0.5 0 0 2];
 
% Set sigma value
sigma = 1e-4;

% compute function value
g = x + sigma* (x'*Q*x)*Q*x;

end
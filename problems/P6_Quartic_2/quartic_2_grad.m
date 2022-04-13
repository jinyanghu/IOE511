function [g] = quartic_2_grad(x)
%QUARTIC_2_GRAD Summary of this function goes here
% Matrix Q
Q = [5 1 0 0.5;
     1 4 0.5 0;
     0 0.5 3 0;
     0.5 0 0 2];
 
% Set sigma value
sigma = 1e4;

% compute function value
g = x + sigma*Q*x;
end


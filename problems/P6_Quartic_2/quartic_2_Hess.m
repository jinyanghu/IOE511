function [H] = quartic_2_Hess(x)
%QUARTIC_2_HESS Summary of this function goes here
%   Detailed explanation goes here
Q = [5 1 0 0.5;
     1 4 0.5 0;
     0 0.5 3 0;
     0.5 0 0 2];
 
% Set sigma value
sigma = 1e4;

% compute function value
H = ones(max(size(x)),1) + 2*sigma*(Q*x)*(x'*Q) + sigma*(x'*Q*x)*Q;
end


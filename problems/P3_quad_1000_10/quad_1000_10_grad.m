function [g] =quad_1000_10_grad(x)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
rng(0);

% Generate random data
q = randn(1000,1);
% MATLAB sprandsym function. Inputs: n, density, reciprocal of the 
% condition number, and kind 
% (see https://www.mathworks.com/help/matlab/ref/sprandsym.html)
Q = sprandsym(1000,0.5,0.1,1);

% compute function value
g = Q*x + q;
end
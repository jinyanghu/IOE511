function [H] = quad_1000_10_Hess(x)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
rng(0);

% Generate random data
q = randn(1000,1);
% MATLAB sprandsym function. Inputs: n, density, reciprocal of the 
% condition number, and kind 
% (see https://www.mathworks.com/help/matlab/ref/sprandsym.html)
Q = sprandsym(1000,0.5,1e-3,1);

% compute function value
H = Q;
end
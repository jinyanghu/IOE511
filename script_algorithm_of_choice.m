% IOE 511/MATH 562, University of Michigan
% Code written by: Albert S. Berahas / Revised by: Jin-Yang Hu

% Script to run code

% close all figures, clear all variables from workspace and clear command
% window
close all
clear all
clc

% set problem (minimal requirement: name of problem)
% Possible problems here, please type string inside ''
% 1. 'Problem7': Rosenbrock_2
% 2. 'Problem8': Rosenbrock_100
problem.name = 'Problem7';

% set problem initial starting point
if problem.name == "Problem7"
    problem.x0 = [-1.2;1];
elseif problem.name == "Problem8"
    problem.x0 = ones(100,1);
    problem.x0(1) = -1.2;
end


% dimension of the problem
problem.n = length(problem.x0);

% final choice of algorithm
% BFGS with backtracking line search
method.name =  'BFGS';
method.options.step_type = 'Backtracking';      


% Optimal choice of parameters for Backtracking line search
% inital_step_size(alpha_hat)
% initial_constant(c) 
% rho(decreasing_factor)
method.options.initial_step_size = 1;
method.options.initial_constant = 1e-8;
method.options.rho = 0.5;

% set options
% 1. norm of the gradient: options.term_tol
% 2. maximum number of iterations: max_iterations
options.term_tol = 1e-6;
options.max_iterations = 1e3;


% run method and return x^* and f^* and (difference between f - f^*, if
% needed)

[x1,f] = optSolver_Hu_Jinyang(problem,method,options);
disp('Done!')
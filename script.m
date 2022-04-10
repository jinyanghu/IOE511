% IOE 511/MATH 562, University of Michigan
% Code written by: Albert S. Berahas / Revised by: Jin-Yang Hu

% Script to run code

% close all figures, clear all variables from workspace and clear command
% window
close all
clear all
clc

% load quadratic data if needed
% load('quadratic10.mat');
% load('quadratic2.mat');

% set problem (minimal requirement: name of problem)
% Possible problems here, please type string inside ''
% 1. 'Quadratic2' / 'Quadratic10'
% 2. 'Rosenbrock'
% 3. 'Function2'
% 4. 'Function3'
problem.name = 'Quadratic2';

if problem.name == "Quadratic2" || problem.name == "Quadratic10"
    if problem.name == "Quadratic2"
        load("quadratic2.mat");
    elseif problem.name == "Quadratic10"
        load("quadratic10.mat");
    end    
    problem.A = A;
    problem.b = b;
    problem.c = c;   
end

if problem.name == "Function2"
    problem.y = [1.5;2.25;2.625];
end

% set problem type: two types of problems 
if problem.name == "Quadratic2" || problem.name == "Quadratic10"
    problem.x0 = x_0; 
elseif problem.name == "Rosenbrock"
    problem.x0 = [1.2;1.2];%[-1.2,1]
elseif problem.name == "Function2"
    problem.x0 = [1;1];
elseif problem.name == "Function3"
    % dimension of the problem, here n = 2, 10, 100, 1000
    % to adjust n to decide dimension of the problem
    n = 1000;
    problem.x0 = zeros(n,1);
    problem.x0(1) = 1;
end


if problem.name == "Rosenbrock"
    problem.x_star = [1;1]; %x_star;
elseif problem.name == "Quadratic2" || problem.name == "Quadratic10"
    problem.x_star = x_star;
end

% dimension of the problem
problem.n = length(problem.x0);

% set method (minimal requirement: name of method)
% possible methods: 
% 1. GradientDescent
% 2. Newton
% 3. BFGS
% 4. LBFGS
% 5. DFP
% 6. TRNewton
% 7. TRSR1
method.name =  'TRNewton';

% Method for each algorithm
% 1. GradientDescent: Backtracking/Constant/Wolfe 
% 2. Newton:  Backtracking/Modified/Wolfe
% 3. BFGS: Backtracking/Wolfe
% 4. LBFGS: Backtracking/Wolfe
method.options.step_type = 'Wolfe';

% step size used by constant GD
method.options.constant_step_size = 1e-3;

% parameters for Backtracking line search
% inital_step_size(alpha_hat)
% initial_constant(c) 
% zro(decreasing_factor)
method.options.initial_step_size = 1;
method.options.initial_constant = 1e-4;
method.options.zro = 0.5;

% parameters for Wolfe line search
method.options.c = 0.5;
method.options.c1 = 1e-8;
method.options.c2 = 1e-4;

%parameters for trust region method
method.options.term_tol_CG = 1e-6;
method.options.region_size = 1;
method.options.tr_c1 = 0.3;
method.options.tr_c2 = 0.8;

% parameters for L-BFGS(number of cuvature pairs)
% n = 2,5,10
% change method.options.m_pairs to decide numbe of curvature pairs used in
% L-BFGS
method.options.m_pairs = 10;

% set options
options.term_tol = 1e-6;
options.max_iterations = 1e3;
options.beta = 1e-6;

% run method and return x^* and f^* and (difference between f - f^*, if
% needed)

[x1,f,f_diff] = optSolver_Hu_Jinyang(problem,method,options);
disp('Done!')



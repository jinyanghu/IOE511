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


% set problem type: two types of problems 

if problem.name == "Problem7"
    problem.x0 = [-1.2;1];
elseif problem.name == "Problem8"
    problem.x0 = ones(100,1);
    problem.x0(1) = -1.2;
end


% dimension of the problem
problem.n = length(problem.x0);

% users can input the method intended to use
% 1. GradientDescent
% 2. GradientDescentW
% 3. Newton
% 4. NewtonW
% 5. TRNewtonCG
% 6. TRSR1CG
% 7. BFGS
% 8. BFGSW
% 9. DFP
% 10. DFPW

method_name = "TRNewtonCG";

if method_name == "GradientDescent"
    method.name = 'GradientDescent';
    method.options.step_type = 'Backtracking';    
elseif method_name == "GradientDescentW"
    method.name =  'GradientDescent';
    method.options.step_type = 'Wolfe';      
elseif method_name == "Newton"
    method.name =  'Newton';
    method.options.step_type = 'Backtracking';      
elseif method_name == "NewtonW"
    method.name =  'Newton';
    method.options.step_type = 'Wolfe';  
elseif method_name == "TRNewtonCG"
    method.name =  'TRNewton';
elseif method_name == "TRSR1CG"
    method.name =  'TRSR1';    
elseif method_name == "BFGS"
    method.name =  'BFGS';
    method.options.step_type = 'Backtracking';      
elseif method_name == "BFGSW"
    method.name = 'BFGS';
    method.options.step_type = 'Wolfe';      
elseif method_name == "DFP"
    method.name = 'DFP';
    method.options.step_type = 'Backtracking';      
elseif method_name == "DFPW"
    method.name = 'DFP';
    method.options.step_type = 'Wolfe';  
else
    error('Method not yet implemented!')
end

% set method (minimal requirement: name of method)
% possible methods: 
% 1. GradientDescent
% 2. Newton
% 3. BFGS
% 4. LBFGS
% 5. DFP
% 6. TRNewton
% 7. TRSR1
%method.name =  'LBFGS';

% Method for each algorithm
% 1. GradientDescent: Backtracking/Constant/Wolfe 
% 2. Newton:  Backtracking/Modified/Wolfe
% 3. BFGS: Backtracking/Wolfe
% 4. LBFGS: Backtracking/Wolfe
% method.options.step_type = 'Backtracking';

% step size used by constant GD
method.options.constant_step_size = 1e-3;

% parameters for Backtracking line search
% inital_step_size(alpha_hat)
% initial_constant(c) 
% rho(decreasing_factor)
method.options.initial_step_size = 1;
method.options.initial_constant = 1e-8;
method.options.rho = 0.5;

% parameters for Wolfe line search
method.options.c = 0.5;
method.options.c1 = 1e-8;
method.options.c2 = 1e-4;

%parameters for trust region method
method.options.term_tol_CG = 1e-6;
methods.options.max_iterations_CG = 1e3;
method.options.region_size = 0.5;
% parameters for TR radius update
method.options.tr_c1 = 0.1;
method.options.tr_c2 = 0.9;


% set options
options.term_tol = 1e-6;
options.max_iterations = 1e3;
options.beta = 1e-6;

% run method and return x^* and f^* and (difference between f - f^*, if
% needed)

[x1,f,f_diff] = optSolver_Hu_Jinyang(problem,method,options);
disp('Done!')
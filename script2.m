% IOE 511/MATH 562, University of Michigan
% Code written by: Albert S. Berahas / Revised by: Jin-Yang Hu

% Script to run code

% close all figures, clear all variables from workspace and clear command
% window
close all
clear all
clc

% load quadratic data if needed
data_p3 = load('problem3.mat');
data_p4 = load('problem3.mat');

% set problem (minimal requirement: name of problem)
% Possible problems here, please type string inside ''
% 1. 'Problem1'
% 2. 'Problem2'
% 3. 'Problem3'
% 4. 'Problem4'
% 5. 'Problem5'
% 6. 'Problem6'
% 7. 'Problem7'
% 8. 'Problem8'
% 9. 'Problem9'
% 10. 'Problem10'
% 11. 'Problem11'
% 12. 'Problem12'


problem.name = 'Problem1';
problems = ["Problem1","Problem2","Problem3","Problem4",...
            "Problem5","Problem6","Problem7","Problem8",...
            "Problem9","Problem10","Problem11","Problem12"];

method_names = ["GradientDescent","GradientDescentW","Newton","NewtonW",...
                "TRNewtonCG","TRSR1CG","BFGS","BFGSW",...
                "DFP","DFPW"];
        


% set problem type: two types of problems 

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
% zro(decreasing_factor)
method.options.initial_step_size = 1;
method.options.initial_constant = 1e-8;
method.options.zro = 0.5;

% parameters for Wolfe line search
method.options.c = 0.5;
method.options.c1 = 1e-8;
method.options.c2 = 1e-4;

%parameters for trust region method
method.options.term_tol_CG = 1e-6;
methods.options.max_iterations_CG = 1e3;
method.options.region_size = 1;
method.options.tr_c1 = 0.1;
method.options.tr_c2 = 0.9;

% parameters for L-BFGS(number of cuvature pairs)
% n = 3,5,8,10
% change method.options.m_pairs to decide numbe of curvature pairs used in
% L-BFGS
method.options.m_pairs = 3;
% strategy to remove curvature pair 
% 1. Random: randomly remove one of curvature pairs in each iteration  
% 2. Min: remove curvature pair with minimum absolute inner product
% 3. Max: remove curvatire pair with maximum absolute inner product
% 4. Oldest: remove oldest curvature
method.strategy = 'Oldest';

% set options
options.term_tol = 1e-6;
options.max_iterations = 1e3;
options.beta = 1e-6;

% run method and return x^* and f^* and (difference between f - f^*, if
% needed)


for i = 1:12
    problem.name = problems(i);
    if problem.name == "Problem1" || problem.name == "Problem2"
        problem.x0 = 20*rand(10,1)-10;
    elseif problem.name == "Problem3" || problem.name == "Problem4"
        problem.x0 = 20*rand(1000,1)-10;
    elseif problem.name == "Problem5" || problem.name == "Problem6"
        problem.x0 = [cos(70) sin(70) cos(70) sin(70)]';
    elseif problem.name == "Problem7"
        problem.x0 = [-1.2;1];
    elseif problem.name == "Problem8"
        problem.x0 = ones(100,1);
        problem.x0(1) = -1.2;
    elseif problem.name == "Problem9"
        problem.x0 = [1;1];
    elseif problem.name == "Problem10"
        n = 10;
        problem.x0 = zeros(n,1);
        problem.x0(1) = 1;
    elseif problem.name == "Problem11"
        n = 1000;
        problem.x0 = zeros(n,1);
        problem.x0(1) = 1;
    elseif problem.name == "Problem12"
        problem.x0 = ones(5,1)*(506.2);
        problem.x0(1) = -506.2;
    end

    if problem.name == "Problem3" || problem.name == "Problem4"
        problem.data = data_p3;
    elseif problem.name == "Problem4"
        problem.data = data_p4;
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
    for j = 1:10
        method_name = method_names(j);

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

        cpu_start_time = cputime;
        [x1,f,f_diff,f_k,g_k] = optSolver_CHL(problem,method,options);
        cpu_time = cputime - cpu_start_time;
        fprintf('Problem %d ,Algorithm: %s Summary Report\n',i,method_names(j));
        fprintf('Iterations: %d\n',length(f_diff));
        fprintf('Function evaluations: %d\n',f_k);
        fprintf('Gradient evaluations: %d\n',g_k);
        fprintf('Cpu seconds: %.4f\n\n',cpu_time);
    end
end
disp('Done!')

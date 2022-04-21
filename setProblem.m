% IOE 511/MATH 562, University of Michigan
% Code written by: Albert S. Berahas / Revised by: Jin-Yang Hu

% Function that specifies the problem. Specifically, a way to compute: 
%    (1) function values; (2) gradients; and, (3) Hessians (if necessary).
%
%           Input: problem (struct), required (problem.name)
%           Output: problem (struct)
%


function [problem] = setProblem(problem)
    
   % calculation of information needed for "quadratic" problem
   % (1) function value: quadratic_f(x)
   % (2) gradient: quadratic_g(x)
   % (3) Hessian: quadtatic_H(x)
   
    function f_value = quadratic_f(x)
        f_value = 1/2 * x' * problem.A * x + problem.b' * x + problem.c;
    end

    function g = quadratic_g(x)
        g = problem.A * x + problem.b;
    end

    function H = quadratic_H(x)
        H = problem.A;
    end
    
    % calculation of information needed for "Rosenbrock"
    % (1) function value: rosen_func(x)
    % (2) gradient: rosen_g(x)
    % (3) Hessian: rosen_H(x)
    function f_value = rosen_func(x)
        f_value = (1 - x(1))^2 + 100*(x(2) - x(1)^2)^2;
    end
    
    function [g] = rosen_g(x)
        g_1 = -2*(1-x(1)) - 400*(x(2)-x(1)^2)*x(1);
        g_2 = 200 * (x(2) - x(1)^2);
        g = [g_1;g_2];
    end
    
    function [H] = rosen_H(x)
        h11 = 2-400*x(2) + 1200 * x(1)^2;
        h22 = 200;
        h12 = -400*x(1);
        H = [h11 h12;h12 h22];
    end

    % calculation of information needed for "Function2"
    % (1) function value: fun2_func(x)
    % (2) gradient: fun2_g(x)
    % (3) Hessian: fun2_H(x)


% check is problem name available
if ~isfield(problem,'name')
    error('Problem name not defined!!!')
end

% set function handles according the the selected problem
switch problem.name

    case 'Problem1'
        problem.compute_f = @quad_10_10_func;
        problem.compute_g = @quad_10_10_grad;
        problem.compute_H = @quad_10_10_Hess;

    case 'Problem2'
        problem.compute_f = @quad_10_1000_func;
        problem.compute_g = @quad_10_1000_grad;
        problem.compute_H = @quad_10_1000_Hess;

    case 'Problem3'
        problem.compute_f = @(x) 1/2*x'*problem.data.Q*x + problem.data.q'*x;
        problem.compute_g = @(x) problem.data.Q*x + problem.data.q;
        problem.compute_H = @(x) problem.data.Q;

    case 'Problem4'
        problem.compute_f = @(x) 1/2*x'*problem.data.Q*x + problem.data.q'*x;
        problem.compute_g = @(x) problem.data.Q*x + problem.data.q;
        problem.compute_H = @(x) problem.data.Q;      

    case 'Problem5'
        problem.compute_f = @quartic_1_func;
        problem.compute_g = @quartic_1_grad;
        problem.compute_H = @quartic_1_Hess;           

    case 'Problem6'
        problem.compute_f = @quartic_2_func;
        problem.compute_g = @quartic_2_grad;
        problem.compute_H = @quartic_2_Hess;         

    case 'Problem7'    
        problem.compute_f = @rosenbrock_2_func;
        problem.compute_g = @rosenbrock_2_grad;
        problem.compute_H = @rosenbrock_2_Hess;            

    case 'Problem8'    
        problem.compute_f = @rosenbrock_100_func;
        problem.compute_g = @rosenbrock_100_grad;
        problem.compute_H = @rosenbrock_100_Hess;  
 
    case 'Problem9'    
        problem.compute_f = @datafit_2_func;
        problem.compute_g = @datafit_2_grad;
        problem.compute_H = @datafit_2_Hess;
        
    case 'Problem10'
        problem.compute_f = @exponential_10_func;
        problem.compute_g = @exponential_10_grad;
        problem.compute_H = @exponential_10_Hess;
    
    case 'Problem11'
        problem.compute_f = @exponential_1000_func;
        problem.compute_g = @exponential_1000_grad;
        problem.compute_H = @exponential_1000_Hess; 
    
    case 'Problem12'
        problem.compute_f = @genhumps_5_func;
        problem.compute_g = @genhumps_5_grad;
        problem.compute_H = @genhumps_5_Hess;         
    
    case 'Quadratic10'
        
        problem.compute_f = @quadratic_f;
        problem.compute_g = @quadratic_g;
        problem.compute_H = @quadratic_H;        
        
    otherwise
        
        error('Problem not defined!!!')
end
end

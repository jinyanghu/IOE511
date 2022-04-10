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

    function f_value = fun2_func(x)
        f_value = 0;
        for i = 1:problem.n
            f_value = f_value + (problem.y(i) - x(1)*(1-x(2)^i))^2;
        end
    end
    
    function [g] = fun2_g(x)
        g1 = 0;
        g2 = 0;
        for i= 1:problem.n
            g1 = g1 + 2*(problem.y(i) - x(1)*(1- x(2)^i))*(-(1-x(2)^i));
            g2 = g2 + 2*(problem.y(i) - x(1)*(1- x(2)^i))*x(1)*i*x(2)^(i-1);
        end
        g = [g1;g2];
    end
    
    function [H] = fun2_H(x)
        h11 = 0;
        h12 = 0;
        h22 = 0;
        for i = 1:problem.n
            h11 = h11+ 2*(problem.y(i) - x(i)*(1-x(2)^i))*(1-x(2)^i)^2;
            h12 = h12 + 2*i*x(2)^(i-1)*(problem.y(i)- x(1)*(1-x(2)^i)-x(1)*((1-x(2)^i)));
            h22 = h22 + 2*x(1)*i*((problem.y(i) - x(1)*(1-x(2)^i))*(i-1)*x(2)^(i-2)+x(1)*i*x(2)^(2*i-2));
        end    
        H = [h11 h12;h12 h22];
    end

    % calculation of information needed for "Function3"
    % (1) function value: fun3_func(x)
    % (2) gradient: fun3_g(x)
    % (3) Hessian: fun3_H(x)

    function f_value = fun3_func(x)
        z1 = (exp(x(1)) - 1)/(exp(x(1)) + 1) + 0.1 * exp(-x(1));
        rest_sum = 0;
        for i = 2:problem.n
            rest_sum = rest_sum + (x(i) - 1)^4;
        end
        f_value = z1 + rest_sum;
    end
    
    function [g] = fun3_g(x)
        g = zeros(problem.n,1);
        g(1) = 2*exp(x(1)) / (exp(x(1)) + 1)^2 - 0.1*exp(-x(1));
        for i = 2:problem.n
            g(i) = 4*(x(i)-1)^3;
        end
    end
    
    function [H] = fun3_H(x)
       H = zeros(problem.n, problem.n);
       H(1,1) = 2*(exp(x(1)) * (exp(x(1)) + 1)^2 -  (exp(x(1)) + 1) * exp(x(1)))  / (exp(x(1)) + 1)^4 + 0.1*exp(-x(1));
       %H(1,1) = 2*(-exp(2*x(1))+exp(x(1)))/(exp(x(1))+1)^3 + 0.1*exp(-x(1));
       for i = 2:problem.n
            H(i,i) = 12*(x(i) - 1)^2;
        end
    end

% check is problem name available
if ~isfield(problem,'name')
    error('Problem name not defined!!!')
end

% set function handles according the the selected problem
switch problem.name

    case 'Function2'
        problem.compute_f = @fun2_func;
        problem.compute_g = @fun2_g;
        problem.compute_H = @fun2_H;

    case 'Function3'
        problem.compute_f = @fun3_func;
        problem.compute_g = @fun3_g;
        problem.compute_H = @fun3_H;

    case 'Rosenbrock'
        
        problem.compute_f = @rosen_func;
        problem.compute_g = @rosen_g;
        problem.compute_H = @rosen_H;
    
    case 'Quadratic2'
        
        problem.compute_f = @quadratic_f;
        problem.compute_g = @quadratic_g;
        problem.compute_H = @quadratic_H;
    
    case 'Quadratic10'
        
        problem.compute_f = @quadratic_f;
        problem.compute_g = @quadratic_g;
        problem.compute_H = @quadratic_H;        
        
    otherwise
        
        error('Problem not defined!!!')
end
end

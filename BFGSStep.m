% IOE 511/MATH 562, University of Michigan
% Code written by: Jin-Yang Hu

% Function that: (1) computes the BFGS step; (2) updates the iterate; and, 
%                (3) computes the function, gradient and hessian at the new iterate
% 
%           Inputs: x, x_old, f, g, g_old, problem, method, options
%           Outputs: x_new, f_new, g_new, d, alpha, skip
%


function [x_new,f_new,g_new,h_new,d,alpha,skip,f_k,g_k] = BFGSStep(x,x_old,f,g,g_old,H,problem,method,options)
% number of evaluations: function/gradient
f_k = 0;
g_k = 0;

% size of x (nx1 vector)
n = size(x,1);

% 0 = not skip, 1 = skip
skip = 0;

% BFGS update for each iteration to decide Hessian approximation
% If it's the first time, set h = identity. Otherwise update h with BFGS method
if class(x_old(1)) == "double" 
    s_k = x - x_old;
    y_k = g - g_old;
    rho = s_k' * y_k;
    % skip update if rho is smaller than the threshold
    if rho <= options.term_tol * norm(s_k) * norm(y_k)
        h = H;
        skip = 1;
    % update if rho is larger than the threshold
    else
        h = (eye(n) - 1/rho * s_k .* y_k') * H * (eye(n) - 1/rho * y_k .* s_k') + 1/rho * (s_k .* s_k');
    end
else
    h = eye(n);
end

% direction
d =  -h * g;

% determine step size
switch method.options.step_type
    case 'Backtracking'
        
        % implementation of Backtracking line search with Armijo condition
        % find alpha satisfies Armijo condition
        alpha = method.options.initial_step_size;
        while problem.compute_f(x + alpha*d) > f + method.options.initial_constant * alpha * g'*d
            alpha = method.options.rho * alpha;
            f_k = f_k + 1;
        end

        % update x,f,g,h
        x_new = x + alpha*d;
        f_new = problem.compute_f(x_new);
        g_new = problem.compute_g(x_new);
        h_new = h;

        % counting times of computing f, g
        f_k = f_k + 1;
        g_k = g_k + 1;

    % implementation of line search with Wolfe condition    
    case 'Wolfe'   
        
        alpha = method.options.initial_step_size;
        alpha_low = 0;
        alpha_high = 1000;

        % find alpha satisfies Wolfe condition
        while true
           x_temp = x + alpha * d;
           f_temp = problem.compute_f(x_temp);
           f_k = f_k + 1;

           % Armijo condition
           threshold = f + method.options.c1 * alpha * g' * d;
           if f_temp <= threshold
               g_temp = problem.compute_g(x_temp);
               g_k = g_k + 1;

               % Curverture condition
               if g_temp' * d >= method.options.c2 * g' * d
                   break;
               end
               alpha_low = alpha;
           else
               alpha_high = alpha;
           end
           alpha = method.options.c*alpha_low + (1-method.options.c)*alpha_high;
        end
        
        % update x,f,g
        x_new = x + alpha*d;
        f_new = problem.compute_f(x_new);
        g_new = problem.compute_g(x_new);
        h_new = h;

        % counting times of computing f, g
        f_k = f_k + 1;
        g_k = g_k + 1;
    otherwise
            error('Method not implemented yet!')
end
end

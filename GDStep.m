% IOE 511/MATH 562, University of Michigan
% Code written by: Albert S. Berahas / Revised by Jin-Yang Hu

% Function that: (1) computes the GD step; (2) updates the iterate; and, 
%                (3) computes the function and gradient at the new iterate
% 
%           Inputs: x, f, g, problem, method, options
%           Outputs: x_new, f_new, g_new, d, alpha
%
function [x_new,f_new,g_new,d,alpha,f_k,gd_k] = GDStep(x,f,g,problem,method,options)

f_k = 0;

% search direction is -g
d = -g;
gd_k = 0;
% determine step size
switch method.options.step_type
    case 'Constant'
        alpha = method.options.constant_step_size;
        x_new = x + alpha*d;
        f_new = problem.compute_f(x_new);
        g_new = problem.compute_g(x_new);
        gd_k = gd_k + 1;
        f_k = f_k + 1;
    case 'Backtracking'
        % implementation of Backtracking line search with Armijo condition        
        alpha = method.options.initial_step_size;
        f_x = problem.compute_f(x);
        while true
           x_temp = x + alpha * d;
           threshold = f_x + method.options.initial_constant * alpha * g' * d;
           if problem.compute_f(x_temp) <= threshold
               f_k = f_k + 1;
               break
           end
           alpha = method.options.zro * alpha;
           f_k = f_k+1;
        end
        
        % fprintf('%.4f, (%.4f,  %.4f) \n',alpha, x(1),x(2));
        x_new = x + alpha*d;
        f_new = problem.compute_f(x_new);
        g_new = problem.compute_g(x_new);
        gd_k = gd_k + 1;
        f_k = f_k+1;
    case 'Wolfe'

        alpha = method.options.initial_step_size;
        alpha_low = 0;
        alpha_high = 1000;
        f_x = problem.compute_f(x);
        g_x = problem.compute_g(x);
        
        while true
           x_temp = x + alpha * d;
           f_temp = problem.compute_f(x_temp);
           f_k = f_k + 1;
           threshold = f_x + method.options.c1 * alpha * g' * d;
           if f_temp <= threshold
               g_temp = problem.compute_g(x_temp);
               gd_k = gd_k + 1;
               if g_temp' * d >= method.options.c2 * g_x' * d
                   break;
               end
               alpha_low = alpha;
           else
               alpha_high = alpha;
           end
           alpha = method.options.c*alpha_low + (1-method.options.c)*alpha_high;
        end
        
        % fprintf('%.4f, (%.4f,  %.4f) \n',alpha, x(1),x(2));
        x_new = x + alpha*d;
        f_new = problem.compute_f(x_new);
        f_k = f_k + 1;
        g_new = problem.compute_g(x_new);
        gd_k = gd_k + 1;
    otherwise
            error('Method not implemented yet!')           
end
end


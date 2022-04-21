% IOE 511/MATH 562, University of Michigan
% Code written by: Jin-Yang Hu

% Function that: (1) computes the BFGS step; (2) updates the iterate; and, 
%                (3) computes the function, gradient and hessian at the new iterate
% 
%           Inputs: x, x_old, f, g, g_old, problem, method, options
%           Outputs: x_new, f_new, g_new, d, alpha, skip
%


function [x_new,f_new,g_new,h_new,d,alpha,skip,f_k,g_k] = DFPStep(x,x_old,f,g,g_old,H,problem,method,options)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% number of evaluations
f_k = 0;
g_k = 0;

skip = 0;
% BFGS update for each iteration to decide Hessian approximation
gd = g;
if class(x_old(1)) == "double" 
    s_k = x - x_old;
    y_k = g - g_old;
    zro = s_k' * y_k;
    if zro <= options.term_tol * norm(s_k) * norm(y_k)
        h = H;
        skip = 1;
    else
        identity = eye(max(size(x)));
        h = H - (H*y_k*y_k'*H)/(y_k'*H*y_k) + s_k*s_k'/(s_k' * y_k);
    end
else
    h = eye(max(size(x)));
end

d =  -h * gd;  

% determine step size
switch method.options.step_type
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
           f_k = f_k + 1;
        end
        
        x_new = x + alpha*d;
        f_new = problem.compute_f(x_new);
        g_new = problem.compute_g(x_new);
        h_new = h;        
        
        f_k = f_k + 1;
        g_k = g_k + 1;
        
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
               g_k = g_k + 1;
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
        g_new = problem.compute_g(x_new);
        h_new = h;
        f_k = f_k + 1;
        g_k = g_k + 1;

    otherwise
            error('Method not implemented yet!')        
                  
        % fprintf('%.4f\n',alpha);        
end
end
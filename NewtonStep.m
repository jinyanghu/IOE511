% IOE 511/MATH 562, University of Michigan
% Code written by: Jin-Yang Hu

% Function that: (1) computes the Newton step; (2) updates the iterate; and, 
%                (3) computes the function, gradient and hessian at the new iterate
% 
%           Inputs: x, f, g, problem, method, options
%           Outputs: x_new, f_new, g_new, d, alpha
%



function [x_new,f_new,g_new,h_new,d,alpha,f_k,g_k] = NewtonStep(x,f,g,H,problem,method,options)

% number of evaluations
f_k = 0;
g_k = 0;

% search direction is -h\g
gd = g;
h = H;
% d = -h\gd;


% determine step size
switch method.options.step_type

    case 'Backtracking'
        
        % modification of Hessian if needed, using Chol decomposition
        d = 0;
        alpha = 1;
        % find smallest Aii
        min_h_element = 1e10;
        h_size = problem.n;
        for i = 1:h_size
            if h(i,i) < min_h_element
                min_h_element = h(i,i);
            end
        end

        identity = eye(h_size);

        if min_h_element > 0
            eta = 0;
        else
            eta = - min_h_element + options.beta;
        end
        
        % make sure that Hessian is PD
        while true
            [R, b] = chol(h+eta * identity);
            if b == 0
                d = -(h+eta * identity)\gd;
                break;
            end
            eta = max(2*eta, options.beta);
        end
        
        % backtracking line search to find alpha for d_k
        alpha = method.options.initial_step_size;
        while problem.compute_f(x + alpha*d) > f + method.options.initial_constant * alpha * g'*d
            alpha = method.options.rho * alpha;
            f_k = f_k + 1;
        end

        x_new = x + alpha*d;
        f_new = problem.compute_f(x_new);
        g_new = problem.compute_g(x_new);
        h_new = problem.compute_H(x_new);
        f_k = f_k + 1;
        g_k = g_k + 1;
    
    case 'Wolfe'
        % modification of Hessian if needed, using Chol decomposition
        d = 0;
        alpha = 1;
        % find smallest Aii
        min_h_element = 1e15;
        h_size = problem.n;
        for i = 1:h_size
            if h(i,i) < min_h_element
                min_h_element = h(i,i);
            end
        end

        identity = eye(h_size);

        if min_h_element > 0
            eta = 0;
        else
            eta = - min_h_element + options.beta;
        end
        
        % make sure that Hessian is PD
        while true
            [R, b] = chol(h+eta * identity);
            if b == 0
                d = -(h+eta * identity)\gd;
                break;
            end
            eta = max(2*eta, options.beta);
        end

        % wolfe line search for step size
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
        h_new = problem.compute_H(x_new);
        f_k = f_k + 1;
        g_k = g_k + 1;

    otherwise
            error('Method not implemented yet!')        
                  
        % fprintf('%.4f\n',alpha);        
end
end



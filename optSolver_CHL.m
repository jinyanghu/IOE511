% IOE 511/MATH 562, University of Michigan
% Code written by: Albert S. Berahas / Revised by: Jin-Yang Hu

% Function that runs a chosen algorithm on a chosen problem
%           Inputs: problem, method, options (structs)
%           Outputs: final iterate (x) and final function value (f),
%           f(x) - f^* at eack iteration
%


function [x,f,f_diff,f_k,gd_k] = optSolver_CHL(problem,method,options)

% set problem, method and options
[problem] = setProblem(problem);
[method] = setMethod(method);
[options] = setOptions(options);

has_x_star = isfield(problem, 'x_star');
%optimal function information
if has_x_star 
    problem.opt_f = problem.compute_f(problem.x_star);
else
    problem.opt_f = 0;
end

% numbe of function evaluations
f_k = 0;
% number of gradient evaluations
gd_k = 0;

% compute initial function/gradient/Hessian
x = problem.x0;
f = problem.compute_f(x);
f_k = f_k + 1;
g = problem.compute_g(x);
gd_k = gd_k+1;

% Hessian initialization for different method
if method.name ~= "BFGS" || method.name~= "LBFGS" || method.name ~= "DFP" ...
   || method.name ~= "TRSR1"     
    [H] = problem.compute_H(x);
else
    [H] = eye(max(size(x)));
end    
norm_g = norm(g,inf);
norm_g0 = norm_g;

% set initial iteration counter
k = 0;

x_old = false;
g_old = false;

% skipped iterations by BFGS
total_skip = 0;

% to store curvature pairs if 'LBFGS' is used
if method.name == "LBFGS"
    curvature_pair = {};
end

% if trust region method is used, parameters needed
if method.name == "TRNewton" || method.name == "TRSR1"
    delta = method.options.region_size;
end

% two termination condition (1) by norm (2) by iterations
while norm_g >= options.term_tol*max(1,norm_g0)  && k < options.max_iterations 
    
    % take step according to a chosen method
    switch method.name
        case 'GradientDescent'
            [x_new,f_new,g_new,d,alpha,if_k,g_k] = GDStep(x,f,g,problem,method,options);
        case 'Newton'
            [x_new,f_new,g_new,h_new,d,alpha,if_k,g_k] = NewtonStep(x,f,g,H,problem,method,options);
        case 'BFGS'
            [x_new,f_new,g_new,h_new,d,alpha,skip,if_k,g_k] = BFGSStep(x,x_old,f,g,g_old,H,problem,method,options);
            total_skip = total_skip + skip;
        case 'DFP'
            [x_new,f_new,g_new,h_new,d,alpha,skip,if_k,g_k] = DFPStep(x,x_old,f,g,g_old,H,problem,method,options);            
        case 'LBFGS'
            [x_new,f_new,g_new,d,alpha,curvature_pair] = ... 
                LBFGSStep(x,x_old,f,g,g_old,curvature_pair,problem,method,options);
                % curvature_pair
        case 'TRNewton'
            [x_new,f_new,g_new,h_new,d,delta_new,if_k,g_k] = NewtonTR(x,f,g,H,delta,problem,method,options);
        case 'TRSR1'
            [x_new,f_new,g_new,h_new,d,delta_new,if_k,g_k] = SR1TR(x,x_old,f,g,g_old,H,delta,problem,method,options);
        otherwise
            
            error('Method not implemented yet!')
            
    end
    % update old and new function values
    x_old = x; f_old = f; g_old = g; norm_g_old = norm_g; 
    x = x_new; f = f_new; g = g_new; norm_g = norm(g,inf);
    gd_k = gd_k + g_k;
    f_k = f_k + if_k;
    % only calculate H if Newton's method is used
    if method.name == "Newton" || method.name == "BFGS" ||  method.name == "DFP"...
       || method.name == "TRNewton" || method.name == "TRSR1" 
        h_old = [H];
        [H] = h_new;
    end
    
    % update delta for TR method
    if method.name == "TRNewton" || method.name == "TRSR1"
        delta_old = delta;
        delta = delta_new;
    end

    % increment iteration counter
    k = k + 1;
    
    % update f(x) - f^* 
    f_diff(k) = problem.compute_f(x_old) - problem.opt_f;
    
    % fprintf('%d, %.4f\n',k,problem.f_diff(k));
end

f_diff(k) = problem.compute_f(x) - problem.opt_f;

% if length of f_diff < k, fill last value in f_diff to k
%if length(f_diff) < options.max_iterations
%    f_end = f_diff(end);
%    f_diff(end+1:options.max_iterations) = f_end;
% end
% total_skip
% fprintf('Total iterations: %d\n',k)
end
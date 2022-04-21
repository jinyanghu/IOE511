% IOE 511/MATH 562, University of Michigan
% Code written by: Jin-Yang Hu

% Function that: (1) computes the Newton step; (2) updates the iterate; and, 
%                (3) computes the function, gradient and hessian at the new iterate
% 
%           Inputs: x, x_old, f, g, g_old,H,curvature_pair, problem, method, options
%           Outputs: x_new, f_new, g_new, d, alpha, curvature_pair
%



function [x_new,f_new,g_new,d,alpha,curvature_pair] = LBFGSStep(x,x_old,f,g,g_old,curvature_pair,problem,method,options)

% BFGS update for each iteration to decide Hessian approximation

%initially start with grgadient
gd = g;
% number of curvature pairs
cv_length = length(curvature_pair);
% start update with identity matrix
h = eye(size(problem.n));

% if there is no curvature pairs, start with gradient
% else do LBFGS algorithm
if class(x_old(1)) ~= "double"
    r = gd;
else
    s_k = x - x_old;
    y_k = g - g_old;
    zro = s_k' * y_k;
    if zro > options.term_tol * norm(s_k) * norm(y_k)
        
        
        if cv_length < method.options.m_pairs
            curvature_pair{cv_length + 1} = {s_k,y_k};
        else
            
           switch method.strategy
                case 'Random'
                    rng('shuffle');
                    remove_idx = randi([1 method.options.m_pairs],1);
                    temp_pair = {};
                    temp_length = 1;
                    for i = 1:cv_length
                        if i == remove_idx
                            continue;
                        else
                            temp_pair{temp_length} = curvature_pair{i};
                            temp_length = temp_length + 1;
                        end
                    end
                    curvature_pair = temp_pair;
                case 'Max'
                    max_idx = 1;
                    max_value = abs(curvature_pair{1}{1}' * curvature_pair{1}{2});
                    for i = 2:cv_length
                        cp_value = abs(curvature_pair{i}{1}' * curvature_pair{i}{2});
                        if cp_value > max_value
                            max_value = cp_value;
                            max_idx = i;
                        end
                    end
                    
                    remove_idx = max_idx;
                    temp_pair = {};
                    temp_length = 1;
                    for i = 1:cv_length
                        if i == remove_idx
                            continue;
                        else
                            temp_pair{temp_length} = curvature_pair{i};
                            temp_length = temp_length + 1;
                        end
                    end
                    curvature_pair = temp_pair;
                case 'Min'
                    min_idx = 1;
                    min_value = abs(curvature_pair{1}{1}' * curvature_pair{1}{2});
                    for i = 2:cv_length
                        cp_value = abs(curvature_pair{i}{1}' * curvature_pair{i}{2});
                        if cp_value < min_value
                            min_value = cp_value;
                            min_idx = i;
                        end
                    end                 

                    remove_idx = min_idx;
                    temp_pair = {};
                    temp_length = 1;
                    for i = 1:cv_length
                        if i == remove_idx
                            continue;
                        else
                            temp_pair{temp_length} = curvature_pair{i};
                            temp_length = temp_length + 1;
                        end
                    end
                    curvature_pair = temp_pair;                    
                
               case 'Oldest'
                    curvature_pair(1) = [];
               otherwise
                   error('Strategy not yet implemented!')
            end
                    
            %curvature_pair(1) = [];
            curvature_pair{cv_length} = {s_k,y_k};
        end
    end    
end

% first loop: start with the most recent pairs
alpha_list = cell(1,cv_length);
zro_list = cell(1,cv_length);
for i = cv_length:-1:1
    si = curvature_pair{i}{1};
    yi = curvature_pair{i}{2};
    zro = 1/ (si' * yi);
    alpha = zro * si' * gd;
    alpha_list{i} = alpha;
    zro_list{i} = zro;
    gd = gd - alpha * yi;
end


% second loop: start with the oldest pairs 
r = h * gd;

for i= 1:cv_length%cv_length:-1:1
    si = curvature_pair{i}{1};
    yi = curvature_pair{i}{2};    
    beta = zro_list{i} * yi' * r;
    r = r + si * (alpha_list{i} - beta);
end

% determine step size
switch method.options.step_type
    case 'Backtracking'
        
        d =  -r;        
        % implementation of Backtracking line search with Armijo condition
        alpha = method.options.initial_step_size;
        f_x = problem.compute_f(x);
        while true
           x_temp = x + alpha * d;
           threshold = f_x + method.options.initial_constant * alpha * g' * d;
           if problem.compute_f(x_temp) <= threshold
               break
           end
           alpha = method.options.zro * alpha;
        end
        
        x_new = x + alpha*d;
        f_new = problem.compute_f(x_new);
        g_new = problem.compute_g(x_new);
        % h_new = problem.compute_H(x_new);          
    
    otherwise
            error('Method not implemented yet!')        
          
        
        % fprintf('%.4f\n',alpha);        
end
end
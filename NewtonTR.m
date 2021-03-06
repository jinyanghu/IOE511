% IOE 511/MATH 562, University of Michigan
% Code written by: Jin-Yang Hu

% Function that: (1) computes the trust region direction; (2) updates the iterate; and, 
%                (3) computes the function, gradient, hessian and size of trust region at 
%                    the new iterate
% 
%           Inputs: x, f, g, H, delta, problem, method, options
%           Outputs: x_new, f_new, g_new, h_new, d, delta_new, f_k, g_k
%

function [x_new,f_new,g_new,h_new,d,delta_new,f_k,g_k] = NewtonTR(x,f,g,H,delta,problem,method,options)
%   Detailed explanation goes here

% number of evaluations: function/gradient
f_k = 0;
g_k = 0;


% implementation of Steihaug CG to solve subproblem
% B: Hessian of objective function at x
B = H;

% initialize parametes
region_size = delta;
z = zeros(problem.n,1);
r = g;
p = -r;

% if norm of r is small, let the direction d = 0
if norm(r,"inf") < method.options.term_tol_CG
    d = z;
end

% CG Steihaug algorithm
while true
    % check if B is non-PD, if so, if so, calculate tao and return || d || = || z + tao*p || = || delta ||
    hessian_p = p'*B*p;
    if hessian_p <= 0
        p_square = norm(p)^2;
        z_square = norm(z)^2;
        tao = (-2*z'*p+sqrt(4*p_square*z_square-4*p_square*(z_square - region_size^2)))/(2*p_square);
        d = z + tao * p;
        break
    end
    
    % update alpha and z
    alpha = (r' * r)/hessian_p;
    z_next = z + alpha*p;

    % check norm of z is larger than the region size , if so, calculate tao and return || d || = || z + tao*p || = || delta ||
    if norm(z_next) >= region_size
        p_square = norm(p)^2;
        z_square = norm(z)^2;       
        tao = (-2*z'*p+sqrt(4*p_square*z_square-4*p_square*(z_square - region_size^2)))/(2*p_square);  
        d = z + tao * p;
        break
    end
    
    % update r
    r_next = r + alpha * B * p;

    % if residual r is samll, return d = z
    if norm(r_next,"inf") <= method.options.term_tol_CG
        d = z_next;
        break
    end

    % update parameters
    beta = (r_next' * r_next)/(r' * r);
    p = -r_next + beta * p;
    r = r_next;
    z = z_next;


end

% trust region method
% evaluation rho to see if the model represents the true cost function well
% if greater than tr_c1, then update
% if greater than tr_c2, then expand trust region (allow to take a larger step)
% else: do not update, then shrink trust region (worse step)
f_new = problem.compute_f(real(x+d));
f_k = f_k + 1;
m_k = f + g' * d + (d'*H*d)/2;
rho = (f - f_new)/(f-m_k);

if rho > method.options.tr_c1
    x_new = x + d;
    g_new = problem.compute_g(x_new);
    h_new = problem.compute_H(x_new);
    g_k = g_k + 1;
    if rho > method.options.tr_c2
        delta_new = 2*delta;
    else
        delta_new = delta;
    end
    x_new = real(x_new);
else
    delta_new = delta/2;
    x_new = x;
    g_new = g;
    f_new = f;
    h_new = H;
end
end
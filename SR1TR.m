function [x_new,f_new,g_new,h_new,d,delta_new] = SR1TR(x,x_old,f,g,g_old,H,delta,problem,method,options)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

B = H;

region_size = delta;
z = zeros(problem.n,1);
r = g;
p = -r;

if norm(r,"inf") < method.options.term_tol_CG
    d = z;
    return 
end

while true
    hessian_p = p'*B*p;
    if hessian_p <= 0
        p_square = norm(p)^2;
        z_square = norm(z)^2;
        tao = (-2*z'*p+sqrt(4*p_square*z_square-4*(z_square - region_size^2)))/(2*p_squre);
        d = z + tao * p;
        break
    end
    
    alpha = (r' * r)/hessian_p;
    z = z + alpha*p;

    if norm(z) >= region_size
        p_square = norm(p)^2;
        z_square = norm(z)^2;
        tao = (-2*z'*p+sqrt(4*p_square*z_square-4*(z_square - region_size^2)))/(2*p_square);        
        d = z + tao * p;
        break
    end
    
    r_next = r + alpha * B * p;
    if norm(r_next,"inf") <= method.options.term_tol_CG
        d = z;
        break
    end

    beta = (r_next' * r_next)/(r' * r);
    p = -r_next + beta * p;
    r = r_next;
end

f_new = problem.compute_f(real(x+d));
m_k = f + g' * d + (d'*H*d)/2;
zro = (f - f_new)/(f-m_k);
if zro > method.options.tr_c1
    x_new = x + d;
    g_new = problem.compute_g(x_new);
    s_k = x_new - x;
    y_k = g_new - g;
    h_new = H + ((s_k - H*y_k)*(s_k - H*y_k)')/((s_k-H*y_k)' * y_k);
    if zro > method.options.tr_c2
        delta_new = 2*delta;
    end
    x_new = real(x_new);    
else
    delta_new = delta/2;
    x_new = x;
    g_new = g;
    h_new = H;
end
end
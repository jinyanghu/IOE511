function [x_new,f_new,g_new,h_new,d,delta_new,f_k,g_k] = SR1TR(x,x_old,f,g,g_old,H,delta,problem,method,options)
% Detailed explanation goes here
% number of evaluations
f_k = 0;
g_k = 0;


if class(x_old(1)) == "double" 
    s_k = x - x_old;
    y_k = g - g_old;
    z = (y_k - H*s_k)' * y_k;
    if abs(z) <= options.term_tol * norm(s_k) * norm((y_k - H*s_k))
        h = H;
    else
        h = H + ((y_k - H*s_k)*(y_k - H*s_k)')/((y_k-H*s_k)' * s_k);
    end
else
    h = H;
end

B = h;

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
        tao1 = (-2*z'*p+sqrt(4*p_square*z_square-4*(z_square - region_size^2)))/(2*p_square);
        g_z = B*z+r;
        tao2 = g_z'*p/hessian_p;
        if abs(tao1) < abs(tao2)
            tao = -tao1;
        else
            tao = tao2;
        end        
        d = z + tao * p;
        break
    end
    
    alpha = (r' * r)/hessian_p;
    z_next = z + alpha*p;

    if norm(z_next) >= region_size
        p_square = norm(p)^2;
        z_square = norm(z)^2;
        tao1 = (-2*z'*p+sqrt(4*p_square*z_square-4*(z_square - region_size^2)))/(2*p_square);        
        g_z = B*z+r;
        tao2 = g_z'*p/hessian_p;

        if abs(tao1) < abs(tao2)
            tao = tao1;
        else
            tao = tao2;
        end       
        d = z + tao * p;
        break
    end
    
    r_next = r + alpha * B * p;
    if norm(r_next,"inf") <= method.options.term_tol_CG
        d = z_next;
        break
    end

    beta = (r_next' * r_next)/(r' * r);
    p = -r_next + beta * p;
    r = r_next;
    z = z_next;
end

f_new = problem.compute_f(real(x+d));
f_k = f_k + 1;
m_k = f + g' * d + (d'*h*d)/2;
rho = (f - f_new)/(f-m_k);

if rho > method.options.tr_c1
    x_new = x + d;
    g_new = problem.compute_g(x_new);
    g_k = g_k + 1;
    if rho > method.options.tr_c2
        delta_new = 2*delta;
    else
        delta_new = delta;
    end
    
    x_new = real(x_new);
    h_new = h;
else
    delta_new = delta/4;
    x_new = x;
    g_new = g;
    h_new = H;
end
end
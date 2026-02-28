function [x,y] = LCP_projection_method(A,b,x0,r,tol,maxiter)
f = @(x) A*x + b;
prox_C = @(z) max(z,0);  %proximal point to R_0^n+
x = projection_method(f,prox_C,x0,r,tol,maxiter);
y = f(x);
y = prox_C(y); %makes solution somewhat nicer
end
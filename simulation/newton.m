function [x,f,dfdx,conv,k] = newton(fun,x0,tol,maxiter)
% Newton's method to find zeros with finite difference approximation of the
% Jacobian
% Remco Leine, INM, University of Stuttgart, 2023
x = x0;
n = length(x0);
k = 0;
conv = 0;
eps = tol/1000;
I = eye(n);
dfdx = zeros(n);
while ~conv && k < maxiter
    f = feval(fun,x);
    if norm(f) < tol
        conv = 1;
    else
        for i = 1:n
            dfdx(:,i) = (feval(fun,x+eps*I(:,i)) - f)/eps;
        end
        k = k + 1;
        x = x - dfdx\f;
    end
end
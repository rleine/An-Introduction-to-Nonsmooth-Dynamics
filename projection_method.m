function x = projection_method(f,proxC,x0,r,tol,maxiter)
% iterative projection method to solve a
% normal cone inclusion problem -f(x) \in N_C(x)
converged = 0;
x = x0;
k = 0;
while ~converged && k<maxiter
    x_old = x;
    x = proxC(x - r*f(x));
    error = norm(x - x_old);
    converged = error < tol;
    k = k + 1;
end
if ~converged, disp('not converged'), end
end
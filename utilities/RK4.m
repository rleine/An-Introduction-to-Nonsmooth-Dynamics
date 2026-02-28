function [t,x] = RK4(odefun,tspan,x0,N)
% Runge-Kutta 4 scheme with fixed stepsize
%
% Input:
%  odefun  right-hand side of the ODE
%  tspan   time span
%  x0      initial condition
%  N       number of time instants
%
% Output:
%  t      time
%  x      state

t0 = tspan(1);
te = tspan(2);
dt = (te-t0)/(N-1);  
t = linspace(t0,te,N);
x = zeros(length(x0),N);
x(:,1) = x0;

for k = 1:(N-1)
    k1 = odefun(t(k),x(:,k));
    k2 = odefun(t(k) + dt/2, x(:,k) + dt/2*k1);
    k3 = odefun(t(k) + dt/2, x(:,k) + dt/2*k2);
    k4 = odefun(t(k) + dt  , x(:,k) + dt*k3);
    x(:,k+1) = x(:,k) + dt/6*(k1 + 2*k2 + 2*k3  + k4);
end

x = x'; 
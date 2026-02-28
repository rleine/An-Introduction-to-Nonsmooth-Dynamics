function f = blockonbelt_arctan(t,x)
% block-on-belt system with arctangent regularization
k = 1;
m = 1;
lambdaN = 10;
mu0 = 0.1;
vdr = 0.2;
delta = 3;
gamma_T = x(2) - vdr;
mu = mu0/(1+delta*abs(gamma_T));
epsilon = 1e3; %regularization parameter
F_T = -mu*lambdaN*2/pi*atan(epsilon*gamma_T);

f = [x(2);
     -k/m*x(1) + F_T/m]; 

      
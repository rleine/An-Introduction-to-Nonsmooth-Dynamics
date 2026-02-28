function f = blockonbelt_LuGre(t,x)
% block-on-belt system with LuGre friction model
k = 1;
m = 1;
lambdaN = 10;
mu0 = 0.1;
vdr = 0.2;
delta = 3;
gamma_T = x(2) - vdr;
mu = mu0/(1+delta*abs(gamma_T));
sigma_0 = 100; %LuGre regularization parameters
sigma_1 = 20;

z = x(3);
zdot = gamma_T -sigma_0*abs(gamma_T)/(mu*lambdaN)*z;
F_T = -(sigma_0*z + sigma_1*zdot);

f = [x(2);
     -k/m*x(1) + F_T/m;
     zdot]; 
      
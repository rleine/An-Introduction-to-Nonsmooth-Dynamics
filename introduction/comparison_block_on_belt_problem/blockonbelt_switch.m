function f = blockonbelt_switch(t,x)
% block-on-belt system with switch model
k = 1;
m = 1;
lambdaN = 10;
mu0 = 0.1;
vdr = 0.2;
delta = 3;
gamma_T = x(2) - vdr;
mu = mu0/(1+delta*abs(gamma_T));
eta = 1e-6; % switch model thickness parameter

if abs(gamma_T)> eta
    F_T = -mu*lambdaN*sign(gamma_T);
    f = [x(2);
          -k/m*x(1) + F_T/m]; 
elseif abs(k*x(1))>mu0*lambdaN
    f = [x(2);
          -k/m*x(1) + mu0*lambdaN/m*sign(k*x(1))]; 
else
    f = [vdr;
         -gamma_T*sqrt(k/m)]; 
end

      
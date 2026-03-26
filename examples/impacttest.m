% impacttest.m , Example 3.5

m = 1; R = 0.1; gravity = 10; mu = 0.05;            % set parameters
eN = 0.8; eT = 0; JS = 1/2*m*R^2;

sys = sys_planar_bouncing_ball(m,JS,R,gravity,mu,eN,eT);  % create system
 
v = 1;                                              % vertical initial speed
Omega_c = 3*(1+eN)*mu*v/(1+eT)/R                    % critical angular vel.

Omega = linspace(-3*Omega_c,3*Omega_c,400);
q = [0; sys.R; 0];                                  % generalized position

for i=1:length(Omega)                               % loop over angular vel.
    umin = [0; -v; Omega(i)];                       % pre-impact gen. vel.
    uplus = impactproblem(sys,0,q,umin,1,1e-8,100); % solve impact probl. 
    phidotplus(i) = uplus(3);                       % post-impact rot. speed
end

plot(Omega,phidotplus)
axis([-6,6,-6,6]), grid
xlabel('$\mathit{\Omega}$'), ylabel('$\dot{\varphi}^+$')
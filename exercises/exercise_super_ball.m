%exercise_super_ball

m = 1;
JS = 0.004;
R = 0.1;
l = 1;
gravity = 10;
mu = 1;
eN = 1;
eT = 1; % Super Ball
sys = sys_super_ball(m,JS,R,l,gravity,mu,eN,eT1);

q0 = [0; 0.5; 0];
u0 = [0.5; 0; 10];

[t,q,u] = moreau_timestepping(sys,0,5,q0,u0,5000,1e-11,1000);

set(figure(1), 'Position', [100, 100, 300, 300]);
plot(q(1,:),q(2,:),'-')
axis([-0.1 0.6 0 0.6]);
xlabel('$x$'), ylabel('$y$')

animate_superball(q, R, 0, 3, 0)
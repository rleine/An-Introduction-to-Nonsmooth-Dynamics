%run_super_ball

m = 1;
JS = 0.004;
R = 0.1;
l = 1;
gravity = 10;
mu = 1;
eN = 1;
eT1 = 1;
eT2 = 0;
sys1 = sys_super_ball(m,JS,R,l,gravity,mu,eN,eT1);
sys2 = sys_super_ball(m,JS,R,l,gravity,mu,eN,eT2);

q0 = [0; 0.9; 0];
u0 = [1; -1; 0];
t0 = 0;
te = 1.2;

%% e_T = 1, ball comes back
[t,q,u] = moreau_timestepping(sys1,t0,te,q0,u0,10000,1e-11,1000);

set(figure(1), 'Position', [100, 100, 300, 300]);
plot(q(1,:),q(2,:),'o-')
axis([-0.1 0.8 0 1]);
title('$e_T = 1$, ball comes back')
xlabel('$x$')
ylabel('$y$')

animate_superball(q, R, 0, 3, 0)


%% e_T = 0, ball does not come back
[t,q,u] = moreau_timestepping(sys2,t0,te,q0,u0,10000,1e-11,1000);

set(figure(3), 'Position', [100, 100, 300, 300]);
plot(q(1,:),q(2,:),'o-')
axis([-0.1 0.8 0 1]);
title('$e_T = 0$, ball does not come back')
xlabel('$x$')
ylabel('$y$')

animate_superball(q, R, 0, 3, 0)
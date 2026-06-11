%run_super_ball

m = 1;
JS = 0.004;
R = 0.1;
l = 1;
gravity = 10;
mu = 1;
eN = 1;
eT1 = 1; % Super Ball
eT2 = 0; % tennis ball
sys1 = sys_super_ball(m,JS,R,l,gravity,mu,eN,eT1);
sys2 = sys_super_ball(m,JS,R,l,gravity,mu,eN,eT2);

q0 = [0; 0.9; 0];
u0 = [1; -1; 0];
t0 = 0;
te = 1.2;

%% e_T = 1, ball comes back
[t_e,q_e,u_e] = eventdriven(sys1,t0,te,q0,u0,1e-11,1e-14,1e-11,1000,1);
set(figure(1), 'Position', [100, 100, 300, 300]);
plot(q_e(1,:),q_e(2,:),'o-')
axis([-0.1 0.8 0 1]);
title('$e_T = 1$, ball comes back')
xlabel('$x$')
ylabel('$y$')

[t,q,u] = moreau_timestepping(sys1,t0,te,q0,u0,10000,1e-11,1000);
animate_superball(q, R, 0, 3, 0)


%% e_T = 0, ball does not come back
[t_e,q_e,u_e] = eventdriven(sys2,t0,te,q0,u0,1e-11,1e-14,1e-11,1000,1);


set(figure(3), 'Position', [100, 100, 300, 300]);
plot(q_e(1,:),q_e(2,:),'o-')
axis([-0.1 0.8 0 1]);
title('$e_T = 0$, ball does not come back')
xlabel('$x$')
ylabel('$y$')

[t,q,u] = moreau_timestepping(sys2,t0,te,q0,u0,10000,1e-11,1000);
animate_superball(q, R, 0, 3, 0)
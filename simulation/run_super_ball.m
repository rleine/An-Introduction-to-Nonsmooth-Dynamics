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

[t,q,u] = eventdriven(sys1,t0,te,q0,u0,1e-11,1e-14,1e-10,1000,1);

fig = figure(1);
set(fig, 'Position', [100, 100, 300, 300]);
plot(q(1,:),q(2,:),'o-')
axis([-0.1 0.8 0 1]);
xlabel('$x$')
ylabel('$y$')


[t,q,u] = eventdriven(sys2,t0,1.4*te,q0,u0,1e-11,1e-14,1e-10,1000,1);

fig = figure(2);
set(fig, 'Position', [100, 100, 300, 300]);
plot(q(1,:),q(2,:),'o-')
axis([-0.1 0.8 0 1]);
xlabel('$x$')
ylabel('$y$')




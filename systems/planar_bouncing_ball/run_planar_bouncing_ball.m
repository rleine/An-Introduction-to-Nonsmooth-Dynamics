% planar bouncing ball 
m = 1; R = 0.1; gravity = 10; mu = 0.05;            % set parameters
eN = 0.8; eT = 0; JS = 0.004;

sys = sys_planar_bouncing_ball(m,JS,R,gravity,mu,eN,eT);  % create system

q0 = [0; 0.1; 0];
u0 = [0; 4.472; -100];

t0 = 0;
te = 7;

[t_e,q_e,u_e] = eventdriven(sys,t0,te,q0,u0,1e-11,1e-14,1e-11,1000,1);
[t_m,q_m,u_m] = moreau_timestepping(sys,t0,te,q0,u0,200,1e-11,1000);
[t_exact,q_exact,u_exact] = moreau_timestepping(sys,t0,te,q0,u0,100000,1e-11,1000);

fig = figure(1);
subplot(2,3,1)
plot(t_e,q_e(1,:),'o-',t_exact,q_exact(1,:))
axis([0 7 0 11]);
xlabel('$t$')
ylabel('$x$')

subplot(2,3,2)
plot(t_e,q_e(2,:),'o-',t_exact,q_exact(2,:))
axis([0 7 -0.5 1.5])
xlabel('$t$')
ylabel('$y$')

subplot(2,3,3)
plot(t_e,q_e(3,:),'o-',t_exact,q_exact(3,:))
axis([0 7 -500 0])
xlabel('$t$')
ylabel('$\varphi$')

subplot(2,3,4)
plot(t_e,u_e(1,:),'o-',t_exact,u_exact(1,:))
axis([0 7 0 3])
xlabel('$t$')
ylabel('$\dot{x}$')

subplot(2,3,5)
plot(t_e,u_e(2,:),'o-',t_exact,u_exact(2,:))
axis([0 7 -10 10])
xlabel('$t$')
ylabel('$\dot{y}$')

subplot(2,3,6)
plot(t_e,u_e(3,:),'o-',t_exact,u_exact(3,:))
axis([0 7 -100 -20])
xlabel('$t$')
ylabel('$\dot{\varphi}$')

% Moreau timestepping method
fig = figure(2);
%set(fig, 'Position', [100, 100, 700, 500]);
subplot(2,3,1)
plot(t_m,q_m(1,:),'o-',t_exact,q_exact(1,:))
axis([0 7 0 11]);
xlabel('$t$')
ylabel('$x$')

subplot(2,3,2)
plot(t_m,q_m(2,:),'o-',t_exact,q_exact(2,:))
axis([0 7 -0.5 1.5])
xlabel('$t$')
ylabel('$y$')

subplot(2,3,3)
plot(t_m,q_m(3,:),'o-',t_exact,q_exact(3,:))
axis([0 7 -500 0])
xlabel('$t$')
ylabel('$\varphi$')

subplot(2,3,4)
plot(t_m,u_m(1,:),'o-',t_exact,u_exact(1,:))
axis([0 7 0 3])
xlabel('$t$')
ylabel('$\dot{x}$')

subplot(2,3,5)
plot(t_m,u_m(2,:),'o-',t_exact,u_exact(2,:))
axis([0 7 -10 10])
xlabel('$t$')
ylabel('$\dot{y}$')

subplot(2,3,6)
plot(t_m,u_m(3,:),'o-',t_exact,u_exact(3,:))
axis([0 7 -100 -20])
xlabel('$t$')
ylabel('$\dot{\varphi}$')

%% Analysis of the Moreau scheme
[t,q,u,PN,PT] = moreau_analysis(sys,t0,te,q0,u0,200,1e-14,1000);

fig = figure(3);
subplot(2,1,1)
plot(t,q(2,:)-sys.R,'o-',t_exact,q_exact(2,:)-sys.R)
axis([0 7 -0.5 1.5])
xlabel('$t$')
ylabel('$g_N$')
grid

subplot(2,1,2)
plot(t,PN,'o-')
axis([0 7 -0.5 10])
xlabel('$t$')
ylabel('$P_N$')
grid

%run_woodpecker

sys = sys_woodpecker;

q0 = [  0;
  -0.000107;
   0.00312];

u0 = [ -0.3505;
        31.81;
        14.44];

[t_e,q_e,u_e] = eventdriven(sys,0,0.146,q0,u0,1e-8,1e-13,1e-10,1000,1);
[t_m,q_m,u_m] = moreau_timestepping(sys,0,0.146,q0,u0,1000,1e-11,1000);
gN_e = sys.getgNall(t_e,q_e);

%% plotting

figure(1); %time histories

subplot(3,1,1)
plot(t_e,q_e(1,:),'o-')
xlabel('$t$')
ylabel('$y$')

subplot(3,1,2)
plot(t_e,q_e(2,:),'o-')
xlabel('$t$')
ylabel('$\varphi_M$')

subplot(3,1,3)
plot(t_e,q_e(3,:),'o-')
xlabel('$t$')
ylabel('$\varphi_S$')


%phase plots
set(figure(2), 'Position', [100, 100, 300, 300]);
plot(q_e(2,:),u_e(2,:),'o-')
axis([-0.15,0.15,-15,35]);
xlabel('$\varphi_M$')
ylabel('$\dot{\varphi}_M$')

set(figure(3), 'Position', [100, 100, 300, 300]);
plot(q_e(3,:),u_e(3,:),'o-')
axis([-0.6,0.2,-15,20]);
xlabel('$\varphi_S$')
ylabel('$\dot{\varphi}_S$')


%% Moreau timestepping

figure(4); %time histories

subplot(3,1,1)
plot(t_m,q_m(1,:),'o-')
xlabel('$t$')
ylabel('$y$')

subplot(3,1,2)
plot(t_m,q_m(2,:),'o-')
xlabel('$t$')
ylabel('$\varphi_M$')

subplot(3,1,3)
plot(t_m,q_m(3,:),'o-')
xlabel('$t$')
ylabel('$\varphi_S$')


%phase plots
set(figure(5), 'Position', [100, 100, 300, 300]);
plot(q_m(2,:),u_m(2,:),'o-')
axis([-0.15,0.15,-15,35]);
xlabel('$\varphi_M$')
ylabel('$\dot{\varphi}_M$')

set(figure(6), 'Position', [100, 100, 300, 300]);
plot(q_m(3,:),u_m(3,:),'o-')
axis([-0.6,0.2,-15,20]);
xlabel('$\varphi_S$')
ylabel('$\dot{\varphi}_S$')


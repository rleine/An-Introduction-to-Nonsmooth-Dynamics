%run_block_on_belt

sys = sys_block_on_belt;

q0 = [0.62; 0];
u0 = [ 0; 0];

%% Quasi-exact solution
[t_exact,q_exact,u_exact] = moreau_timestepping(sys,0,30,q0,u0,100000,1e-8,1000);

%% event-driven
[t_e,q_e,u_e] = eventdriven(sys,0,30,q0,u0,1e-8,1e-8,1e-13,1000,1);

%phase plot
set(figure(1), 'Position', [100, 100, 300, 300]);
plot(q_e(1,:),u_e(1,:),'o-')
axis([-0.6 1.2 -1 0.4]);
xlabel('$x$'); ylabel('$\dot{x}$'); 

%time history
set(figure(2), 'Position', [100, 100, 300, 300]);
plot(t_e,q_e(1,:),'o-')
xlabel('$t$');ylabel('$x$')


%% Moreau timestepping scheme
q0 = [0.62; 0];
[t_m,q_m,u_m] = moreau_timestepping(sys,0,30,q0,u0,8000,1e-8,1000);

%phase plot
set(figure(3), 'Position', [100, 100, 300, 300]);
plot(q_exact(1,:),u_exact(1,:),'k',q_e(1,:),u_e(1,:),'bo-',q_m(1,:),u_m(1,:),'ro-')
axis([-0.6 1.2 -1 0.4]);
xlabel('$x$'); ylabel('$\dot{x}$')

%time history
set(figure(4), 'Position', [100, 100, 300, 300]);
plot(t_exact,q_exact(1,:),'k',t_e,q_e(1,:),'bo-',t_m,q_m(1,:),'ro-')
xlabel('$t$'); ylabel('$x$')


%% event-driven
tic
for i = 1:100
    [t_e,q_e,u_e] = eventdriven(sys,0,30,q0,u0,1e-8,1e-8,1e-13,1000,0);
end
toc

%% Moreau
tic
for i = 1:100
    [t_m,q_m,u_m] = moreau_timestepping(sys,0,30,q0,u0,8000,1e-8,1000);
end
toc



%run_rockingblock
% Exercise 3.2

% Exercise 3.2a
sys = sys_rocking_block(0.0125,0.05);

phi0 = 8/180*pi;
q0 = [0;sys.a*sin(phi0)+sys.b*cos(phi0);phi0];       
u0 = [0;0;0];

[t,q,u] = moreau_timestepping(sys,0,2,q0,u0,4000,1e-13,3000);
gN = sys.getgNall(t,q);

fig = figure(1)
set(fig, 'Position', [100, 100, 300, 300]);
plot(t,gN);
xlabel('$t$'), ylabel('$g_N$')
legend('$g_{N1}$','$g_{N2}$')

% Exercise 3.2b
sys = sys_rocking_block(0.05,0.0125);

phi0 = 30/180*pi;
q0 = [0;sys.a*sin(phi0)+sys.b*cos(phi0);phi0];       
u0 = [0;0;0];

[t,q,u] = moreau_timestepping(sys,0,2,q0,u0,4000,1e-13,3000);
gN = sys.getgNall(t,q);

fig = figure(2)
set(fig, 'Position', [100, 100, 300, 300]);
plot(t,gN);
xlabel('$t$'), ylabel('$g_N$')
legend('$g_{N1}$','$g_{N2}$')


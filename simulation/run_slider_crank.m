% slider crank system 
mu = 0.01;            % set parameters
eN = 0.4; 
eT = 0;

sys = sys_slider_crank(mu,eN,eT);  % create system

q0 = [0; 0; 0.017];
u0 = [150; -75; 0];

t0 = 0;
te = 0.2;

%[t_e,q_e,u_e] = eventdriven(sys,t0,te,q0,u0,1e-13,1e-14,1e-13,1000,1);
[t_m,q_m,u_m] = moreau_timestepping(sys,t0,te,q0,u0,10000,1e-11,1000);
[t_exact,q_exact,u_exact] = moreau_timestepping(sys,t0,te,q0,u0,200000,1e-11,1000);

gN_m = sys.getgNall(t_m,q_m);
%gN_e = sys.getgNall(t_e,q_e);
gN_exact = sys.getgNall(t_exact,q_exact);
fig = figure(1)
set(fig, 'Position', [100, 100, 300, 300]);
plot(t_m,gN_m(3,:),'or',t_exact,gN_exact(3,:),'k')
xlabel('$t$')
ylabel('$g_{N3}$')
grid

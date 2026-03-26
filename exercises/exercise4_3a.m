
a = 0.0125; 
b = 0.05;
m = 1;
J = 1/12*m*(4*a^2+4*b^2);
sys = sys_rocking_block(a,b,m,J,9.81,0.3,1e-1);

phi0 = 8/180*pi;
q0 = [0;sys.a*sin(phi0)+sys.b*cos(phi0);phi0];       
u0 = [0;0;0];

[t,q,u] = moreau_timestepping(sys,0,2,q0,u0,4000,1e-13,3000);
gN = sys.getgNall(t,q);

figure(1)
plot(t,gN);
xlabel('$t$'), ylabel('$g_N$')
legend('$g_{N1}$','$g_{N2}$')
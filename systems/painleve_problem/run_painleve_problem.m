%run_painleve_problem
mu = 5/3;
sys = sys_painleve_problem(mu);

q0 = [0;0.515;31/180*pi];       
u0 = [30;0;0];

[t,q,u] = moreau_timestepping(sys,0,1.5,q0,u0,4000,1e-13,10000);

gN = sys.getgNall(t,q);

[gammaN,gammaT] = deal(zeros(length(t),1));

A = @(varphi) 1/sys.m + (sys.s^2)/sys.JS*cos(varphi).^2.*(1 - mu*tan(varphi));
b = @(varphi,varphidot) -sys.gravity + sys.s*sin(varphi).*varphidot.^2; 

for i=1:length(t)
 [WN,WT] = sys.Wchi(t(i),q(:,i),sys.I);
 gammaN(i) = WN'*u(:,i);
 gammaT(i) = WT'*u(:,i);
end


fig = figure(1);
subplot(2,3,1)
plot(t,u(1,:),'o-');
xlabel('$t$'), ylabel('$\dot{x}$'), grid

subplot(2,3,2)
plot(t,u(2,:),'o-');
xlabel('$t$'), ylabel('$\dot{y}$'), grid


subplot(2,3,3);
plot(t,u(3,:),'o-');
xlabel('$t$'), ylabel('$\dot{\varphi}$'), grid

subplot(2,3,4)
plot(t,gN,'o-');
xlabel('$t$'), ylabel('$g_N$'), grid

subplot(2,3,5)
plot(t,gammaN,'o-');
xlabel('$t$'), ylabel('$\gamma_N$'), grid

subplot(2,3,6)
plot(t,gammaT,'o-');
xlabel('$t$'), ylabel('$\gamma_T$'), grid

figure(2)
plot(t,A(q(3,:)),'o-',t,b(q(3,:),u(3,:)),'o-')
xlabel('$t$')
legend('A','b')
grid

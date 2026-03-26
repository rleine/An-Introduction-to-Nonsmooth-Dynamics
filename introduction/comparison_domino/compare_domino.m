% compare_domino.m
%
% Data from:
% Lu, Third & Müller, Effect of particle shape on domino wave propagation: 
% a perspective from 3D, anisotropic discrete element simulations
% DOI 10.1007/s10035-013-0472-2

run_domino
v = calc_speed(t,phi,n,d);

%% Simulate dominoes using Hunt-Crossley and smoothing
kHC = 5e4;
cHC = 2;
epsilon = 1000;
x0 = [q0;u0];


disp('Simulation of dominoes using Hunt-Crossley contact model and smooth friction RK4')
tic
[t_RK4,x_RK4] = RK4(@(t,x) ode_HC(t,x,sys,kHC,cHC,epsilon),[0 t1],x0,1500);
toc
q_RK4 = x_RK4(:,1:sys.ndof)';
phi_RK4 = q_RK4(3:3:end,:);
v_RK4 = calc_speed(t_RK4,phi_RK4,n,d);

disp('Simulation of dominoes using Hunt-Crossley contact model and smooth friction ode45')
tic
opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t_ode45,x_ode45] = ode45(@(t,x) ode_HC(t,x,sys,kHC,cHC,epsilon),[0 t1],x0,opts);
toc
q_ode45 = x_ode45(:,1:sys.ndof)';
phi_ode45 = q_ode45(3:3:end,:);
v_ode45 = calc_speed(t_ode45,phi_ode45,n,d);

%% post-processing and compare
fig = openfig('lambda=1.8t.fig','invisible');
dataObjs = findobj(fig,'-property','YData');
k_exp = dataObjs(1).XData';
v_exp = dataObjs(1).YData';
close all

figure(1)
plot(1:n,v,'*-r',1:n,v_RK4,'v-b',1:n,v_ode45,'^-m',k_exp,v_exp,'o-k')
xlabel('domino sequence number')
ylabel('intrinsic collision speed [m/s]')
legend('Moreau','Hunt-Crossley RK4','Hunt-Crossley ode45','Measurement data')


fig = openfig('inclined_angle_L=1.8h.fig','invisible');
dataObjs = findobj(fig,'-property','YData');
t_exp = dataObjs(1).XData';
phi_exp = dataObjs(1).YData';
close(fig);

figure(2)
num = 20;
offset = @(phi) find(abs(phi(num, :)) > 1e-3, 1, 'first');
plot(t-t(offset(phi)),phi(num,:)*180/pi+90,'r',...
     t_RK4-t_RK4(offset(phi_RK4)),phi_RK4(num,:)*180/pi+90,'b',...
     t_ode45-t_ode45(offset(phi_ode45)),phi_ode45(num,:)*180/pi+90,'m',...
     t_exp, phi_exp,'o-k')
axis([-0.001 0.22 10 100])
xlabel('Shifted time [s]')
ylabel('Domino inclined angle (deg)')
legend('Moreau','Hunt-Crossley RK4','Hunt-Crossley ode45','Measurement data')

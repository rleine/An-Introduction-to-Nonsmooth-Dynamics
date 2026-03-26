% Comparison of Hunt-Crossley versus Nonsmooth Dynamics approach
% on the rocking block problem
%
% Experimental data from
%
% Shaking Table Seismic Experimental Investigation of
% Lightweight Rigid Bodies
% 
% Giuseppe Cocuzza Avellino, Francesco Cannizzaro, Nicola Impollonia
%
% https://doi.org/10.3390/buildings12070915
% https://zenodo.org/records/6670588#.YrB6FHZBw7c

%% Read dataset

[num, txt, raw] = xlsread('Dataset.xlsx',2);
data = cell2mat(raw(4:end, :));

block = 'B2';

switch block
    case 'A2'
        m = 0.435;
        a = 0.0285;
        b = 0.1423;
        alpha= 0.70; % correction factor
        t_end = 10;
   case 'B1'
        m = 0.995;
        a = 0.0385;
        b = 0.15;
        alpha= 0.739; % correction factor
        t_end = 5;
   case 'B2'
        m = 0.3;
        a = 0.0282;
        b = 0.1065;
        alpha= 0.909; % correction factor
        t_end = 5;
   case 'C2'
        m = 0.201;
        a = 0.0283;
        b = 0.0710;
        alpha= 0.879; % correction factor
        t_end = 2;
end

index = (double(upper(block(1))) - double('A'))*3 + str2double(block(2));
t_e = data(:,index*2-1);
theta = data(:,index*2);
g = 9.801; % constant of gravity in Catania
mu = 0.5;
J = 1/12*m*(4*a^2+4*b^2);

theta_dot = gradient(theta,t_e);
phi0 = theta(1);
phidot0 = theta_dot(1);
ydot0 = -(-a*cos(phi0)*sign(phi0)+b*sin(phi0))*phidot0;
xdot0 = -(a*sin(phi0)*sign(phi0)+b*cos(phi0))*phidot0;

q0 = [0;a*sin(abs(phi0))+b*cos(abs(phi0));phi0];       
u0 = [xdot0;ydot0;phidot0];


%% Simulations using Moreau timestepping scheme; test block B2
sys = sys_rocking_block(a,b,m,J,g,mu,1e-2);
[t_m,q_m,u_m] = moreau_timestepping(sys,0,t_end,q0,u0,30000,1e-9,1e4);

ydot0 = -(-alpha*a*cos(phi0)*sign(phi0)+b*sin(phi0))*phidot0;
xdot0 = -(alpha*a*sin(phi0)*sign(phi0)+b*cos(phi0))*phidot0;
q0 = [0;alpha*a*sin(abs(phi0))+b*cos(abs(phi0));phi0]; 
u0 = [xdot0;ydot0;phidot0];
sys = sys_rocking_block(alpha*a,b,m,J,g,mu,1e-2);
[t_mc,q_mc,u_mc] = moreau_timestepping(sys,0,t_end,q0,u0,30000,1e-9,1e4);

set(figure(1), 'Position', [100, 100, 700, 200]);
plot(t_e,theta,'k',t_m,q_m(3,:),'b--',t_mc,q_mc(3,:),'r-.')
xlabel('$t$'), ylabel('$\varphi$')
legend('measurement','Moreau, theoretical','Moreau, corrected')

set(figure(2), 'Position', [100, 100, 700, 200]);
plot(t_e,theta_dot,'k',t_m,u_m(3,:),'b--',t_mc,u_mc(3,:),'r-.')
xlabel('$t$'), ylabel('$\dot{\varphi}$')
legend('measurement','Moreau, theoretical','Moreau, corrected')


%% Simulations using Hunt-Crossley and regularized friction
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);

k = 1e7; % Hunt-Crossley contact stiffness
epsilon = 1e4;

sys_HC = @(t,x,c) rocking_block_HC(t,x,alpha*a,b,m,J,g,mu,k,c,epsilon);

[t1,x1] = ode45(@(t,x) sys_HC(t,x,10),[0 t_end],[q0;u0],opts);
[t2,x2] = ode45(@(t,x) sys_HC(t,x,40),[0 t_end],[q0;u0],opts);
[t3,x3] = ode45(@(t,x) sys_HC(t,x,100),[0 t_end],[q0;u0],opts);

set(figure(3), 'Position', [100, 100, 700, 200]);
plot(t_e,theta,'k',t1,x1(:,3),'r:',t2,x2(:,3),'b--',t3,x3(:,3),'m-.')
xlabel('$t$'), ylabel('$\varphi$')
labels = {'measurement',...
       'Hunt-Crossley $c_\mathrm{HC} = 10$', ...
       'Hunt-Crossley $c_\mathrm{HC} = 40$', ...
       'Hunt-Crossley $c_\mathrm{HC} = 100$'};
legend(labels)

set(figure(4), 'Position', [100, 100, 700, 200]);
plot(t_e,theta_dot,'k',t1,x1(:,6),'r:',t2,x2(:,6),'b--',t3,x3(:,6),'m-.')
xlabel('$t$'), ylabel('$\dot{\varphi}$')
legend(labels)
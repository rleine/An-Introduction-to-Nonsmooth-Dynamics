% run_blockonbelt_comparison.m
% comparison of time integration methods
% for the block-on-belt system with stick-slip limit cycle

ax = [-0.7 1.2 -0.9 0.3];

% quasi-exact solution using Moreau's timestepping scheme
sys = sys_block_on_belt;
q0 = [1.135; 0];
u0 = [ 0; 0];
[t_ex,q_ex,u_ex] = moreau_timestepping(sys,0,12.3,q0,u0,100000,1e-8,1000);
x_ex= [q_ex(1,:); u_ex(1,:)]';

% regularization using arctangent function
tols = odeset('RelTol',1e-5,'AbsTol',1e-5);
[t_atan,x_atan] = ode45('blockonbelt_arctan',[0 11.3],[1.067; 0],tols);
figure(1)
plot(x_atan(:,1),x_atan(:,2),'o-',x_ex(:,1),x_ex(:,2),'k')
axis(ax); xlabel('$x$'); ylabel('$\dot{x}$')
disp(['arctangent regularization, time steps = ',num2str(length(t_atan))])

% LuGre model
x0 = [1.043656607367126   0.081148750297240  -0.007831711872955]';
tols = odeset('RelTol',1e-5,'AbsTol',1e-5);
[t_LuGre,x_LuGre] = ode45('blockonbelt_LuGre',[0 11.05],x0,tols);
figure(2)
plot(x_LuGre(:,1),x_LuGre(:,2),'o-',x_ex(:,1),x_ex(:,2),'k')
axis(ax); xlabel('$x$'); ylabel('$\dot{x}$')
disp(['LuGre method, time steps = ',num2str(length(t_LuGre))])

% switch model
x0 = [1.135   0]';
tols = odeset('RelTol',1e-6,'AbsTol',1e-6);
[t_switch,x_switch] = ode45('blockonbelt_switch',[0 11.96],x0,tols);
figure(3)
plot(x_switch(:,1),x_switch(:,2),'o-',x_ex(:,1),x_ex(:,2),'k')
axis(ax); xlabel('$x$'); ylabel('$\dot{x}$')
disp(['Switch model, time steps = ',num2str(length(t_switch))])
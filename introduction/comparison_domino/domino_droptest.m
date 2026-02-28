% Drop test

% Data from:
% Lu, Third & Müller, Effect of particle shape on domino wave propagation: a
% perspective from 3D, anisotropic discrete element simulations
% DOI 10.1007/s10035-013-0472-2

n = 1;           % number of dominos
a = 0.5*7.7e-3;   % half-thickness
b = 0.5*43.2e-3;  % half-height
w = 21.9e-3;      % width
rho = 631;        % density


g = 9.80665;      % acceleration of gravity in Zurich 
muG = 0.90;       % friction coefficient with the ground
muD = 0.11;       % friction coefficient between dominos
eG = 0.00;        % restitution coefficient with the ground
eD = 0.54;        % restitution coefficient between dominos

m = rho*2*a*2*b*w;
J = 1/12*m*(4*a^2+4*b^2);

%% Simulation of droptest using Moreau timestepping scheme

sys = sys_domino(n,a,b,m,J,g,muG,muD,eD,eD,1e-5);


q0 = [0;3*b;-pi+atan(a/b)];
u0 = zeros(3,1);

disp('Moreau timesteping simulation of domino drop test')
tic
[t,q,u] = moreau_timestepping(sys,0,0.2,q0,u0,10000,1e-6,1e4);
toc

%% Simulate of droptest using Hunt-Crossley and smoothing
kHC = 5e5;
cHC = 2;
epsilon = 1000;
x0 = [q0;u0];
opts = odeset('RelTol',1e-6,'AbsTol',1e-6);

disp('Simulation of dominoes using Hunt-Crossley contact model and smooth friction')
tic
%[t_p,x_p] = ode45(@(t,x) ode_HC(t,x,sys,kHC,cHC,epsilon),[0 1],x0,opts);
[t_p,x_p] = RK4(@(t,x) ode_HC(t,x,sys,kHC,cHC,epsilon),[0 0.2],x0,1000);
q_p = x_p(:,1:sys.ndof)';
u_p = x_p(:,sys.ndof+1:end)';
%animate_dominos(q_p, n, a, b,0,1);

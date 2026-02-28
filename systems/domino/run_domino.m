% run_domino.m
%
% Data from:
% Lu, Third & Müller, Effect of particle shape on domino wave propagation: 
% a perspective from 3D, anisotropic discrete element simulations
% DOI 10.1007/s10035-013-0472-2

t1 = 1;

n = 30;           % number of dominos
a = 0.5*7.7e-3;   % half-thickness
b = 0.5*43.2e-3;  % half-height
w = 21.9e-3;      % width
rho = 631;        % density
d = 1.8*2*a;      % initial width between dominoes
g = 9.80665;      % acceleration of gravity in Zurich 
muG = 0.90;       % friction coefficient with the ground
muD = 0.11;       % friction coefficient between dominoes
eG = 0.00;        % restitution coefficient with the ground
eD = 0.54;        % restitution coefficient between dominoes

m = rho*2*a*2*b*w;
J = 1/12*m*(4*a^2+4*b^2);

%% Simulation using Moreau timestepping scheme

sys = sys_domino(n,a,b,m,J,g,muG,muD,eG,eD,1e-5);

% initial condition
phi10 = 0;
phidot10 = -12;
ydot10 = -(a*cos(phi10)+b*sin(phi10))*phidot10;
xdot10 = -(a*sin(phi10)+b*cos(phi10))*phidot10;
q10 = [0;-a*sin(phi10)+b*cos(phi10);phi10]; 
u10 = [xdot10;ydot10;phidot10];

q0 = reshape([(d+2*a)*((1:n).'-1), b*ones(n,1), zeros(n,1)].', [], 1);
q0(1:3) = q10;
u0 = zeros(3*n,1);
u0(1:3) = u10;

disp('Moreau timestepping simulation of dominoes')
tic
[t,q,u] = moreau_timestepping(sys,0,t1,q0,u0,4000,1e-6,1e4);
toc

phi = q(3:3:3*n,:);
plot(t,phi)


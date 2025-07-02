% Simulation and animation of the superball 
% Remco Leine, May 2024
clc
clear all
close all
global R l

q0 = [0; 0.9; 0];
u0 = [1; -1; 0];
t0 = 0;
te = 1.2;
Ni = 500;
tol = 1e-6;
maxiter = 1000;
[t,q,u] = timestepping('superball_sys',t0,te,q0,u0,Ni,tol,maxiter);



%% animation with rotation
figure
hold on
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Initializing
  % Initial Ball points
  th=0:.1:2*pi;
  xx=R*cos(th)+q0(1);
  yy=R*sin(th)+q0(2);
  % Initial Ball line points
  x1l1=1*R*cos(q0(3))+q0(1);
  x1l2=1*R*cos(q0(3)+pi)+q0(1);
  y1l1=1*R*sin(q0(3))+q0(2);
  y1l2=1*R*sin(q0(3)+pi)+q0(2);
  
  x2l1=1*R*cos(q0(3)+pi/2)+q0(1);
  x2l2=1*R*cos(q0(3)+pi+pi/2)+q0(1);
  y2l1=1*R*sin(q0(3)+pi/2)+q0(2);
  y2l2=1*R*sin(q0(3)+pi+pi/2)+q0(2);
  
  % Drawing the floor
  floormin=-0.5;
  floorwidth=2;
  floorheight=floorwidth;
  heightmin=-1;
  axis([floormin floormin+floorwidth heightmin heightmin+floorheight]);
  plot(floormin:.1:floormin+floorwidth,zeros(length(floormin:.1:floormin+floorwidth)),'k','LineWidth',2); 


    % Drawing the ceiling
  ceilingmin= 0.3;
  ceilingwidth=1;
  ceilingheight=ceilingwidth;
  heightmin=-1;
  
  plot(ceilingmin:.1:ceilingmin+ceilingwidth,l+zeros(length(ceilingmin:.1:ceilingmin+ceilingwidth)),'k','LineWidth',2);

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Initial drawing of ball
  Ballhndl1=plot([x1l1 x1l2],[y1l1 y1l2],'red');
  Ballhndl2=plot([x2l1 x2l2],[y2l1 y2l2],'red');
  Ballhndl=plot(xx,yy,'blue');%,'EraseMode','background');
%   set(gca,'UserData',[Ballhndl;Ballhndl1;Ballhndl2]);
  
  
  % animation
%         set(0,'currentfigure',Ball)
%       Ballhndl=get(gca,'UserData');
      Ballhndl=[Ballhndl;Ballhndl1;Ballhndl2];
  
  for i=1:Ni    
  %Ball Circle points
  th=0:.1:2*pi;
  x1=q(1,i)+R*cos(th);
  y1=q(2,i)+R*sin(th);
  
     
  %%%%%%%%%%%%%% Ball line points     
  x1l1=1*R*cos(q(3,i))+q(1,i);
  x1l2=1*R*cos(q(3,i)+pi)+q(1,i);
  y1l1=1*R*sin(q(3,i))+q(2,i);
  y1l2=1*R*sin(q(3,i)+pi)+q(2,i);
  
  x2l1=1*R*cos(q(3,i)+pi/2)+q(1,i);
  x2l2=1*R*cos(q(3,i)+pi+pi/2)+q(1,i);
  y2l1=1*R*sin(q(3,i)+pi/2)+q(2,i);
  y2l2=1*R*sin(q(3,i)+pi+pi/2)+q(2,i);
  
  % Draw Lines
   set(Ballhndl(1),'XData',x1); 
   set(Ballhndl(1),'YData',y1);
   set(Ballhndl(2),'XData',[x1l1 x1l2]);
   set(Ballhndl(2),'YData',[y1l1 y1l2]);
   set(Ballhndl(3),'XData',[x2l1 x2l2]);
   set(Ballhndl(3),'YData',[y2l1 y2l2]);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  axis equal
%    drawnow;
    pause(0.0001);
  end
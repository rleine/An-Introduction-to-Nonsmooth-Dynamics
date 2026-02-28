function xdot = rocking_block_HC(t,x,a,b,m,J,g,mu,k,c,epsilon)
% planar rocking block with Hunt-Crossley and smoothed friction
q = x(1:3);
qdot = x(4:end);
M = diag([m,m,J]);
h = [ 0;
      -m*g;
      0];
y = q(2);
phi = q(3);
     
WN(:,1) = [ 0; 1;-a*cos(phi)+b*sin(phi)];
WN(:,2) = [ 0; 1; a*cos(phi)+b*sin(phi)];
WT(:,1) = [ 1; 0; a*sin(phi)+b*cos(phi)];
WT(:,2) = [ 1; 0;-a*sin(phi)+b*cos(phi)];
gN1 = y-a*sin(phi)-b*cos(phi);
gN2 = y+a*sin(phi)-b*cos(phi);
gN =[gN1;gN2];
gammaN = WN'*qdot;
gammaT = WT'*qdot;
delta = max(-gN,0);
deltadot = -gammaN;
F_N = k*delta.^1.5.*(1 + c.*deltadot); %Hunt-Crossley model
F_T = -mu.*F_N.*2./pi.*atan(epsilon*gammaT);
xdot = [qdot;
        M\(h + WN*F_N + WT*F_T)];
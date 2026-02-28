function xdot = ode_HC(t,x,sys,kHC,cHC,epsilon)
ndof = length(x)/2;
q = x(1:ndof);
u = x(ndof+1:end);
M = sys.M(t,q);
h = sys.h(t,q,u);
[WN,WT,chiN,chiT] = sys.Wchi(t,q,sys.I);
gN = sys.contactdistance(t,q,sys.I);
gammaN = WN'*u + chiN;
gammaT = WT'*u + chiT;
delta = max(-gN,0);
deltadot = -gammaN;
FN = kHC*delta.^1.5.*(1 + cHC.*deltadot);
FT = -sys.mu.*FN.*2/pi.*atan(epsilon*gammaT);
xdot = [u;
        M\(h + WN*FN + WT*FT)]; 
function [IP,IPS,IPLplus,IPLmin,udot] = contactaccelerationproblem(sys,...
    t,q,u,IN,IS,ILplus,ILmin,tol,maxiter)
M = sys.M(t,q);
h = sys.h(t,q,u);

% selection matrices 
S = eye(length(sys.I));
SLplus = S(ILplus,IN);  % forward slipping contacts
SLmin = S(ILmin,IN);    % forward slipping contacts
SS = S(IS,IN);          % sticking contacts

[WN,WT] = sys.Wchi(t,q,IN); % generalized force directions
[nuN,nuT] = sys.nu(t,q,u,IN);
WS = WT*SS';
WLplus = WT*SLplus';
WLmin = WT*SLmin';
WL = [WLplus WLmin];
nuS = SS*nuT;

mu = sys.contactpar(IN); % friction coefficients
QL = [- SLplus*diag(mu); 
        SLmin*diag(mu)]; 

% Delassus submatrices
GNN = WN'*(M\WN); GNL = WN'*(M\WL); GNS = WN'*(M\WS);
GSN = WS'*(M\WN); GSL = WS'*(M\WL); GSS = WS'*(M\WS);

A = [GNN+GNL*QL GNS;  % use gammaNdot = gammaSdot = 0 as first guess
     GSN+GSL*QL GSS]; % by solving the system A*x + b = 0

b = [WN'*(M\h) + nuN;
     WS'*(M\h) + nuS];

x = -A\b;

lambdaN = x(1:length(IN));
lambdaS = x(length(IN)+1:end,1);

converged = 0;
k = 0;
while ~converged && k<maxiter
    lambdaN_old = lambdaN;
    lambdaS_old = lambdaS;
    gammaNdot = (GNN+GNL*QL)*lambdaN + GNS*lambdaS + WN'*(M\h) + nuN;
    gammaSdot = (GSN+GSL*QL)*lambdaN + GSS*lambdaS + WS'*(M\h) + nuS;
    lambdaN = -prox_CN(-lambdaN + sys.rN*gammaNdot);
    lambdaS = -prox_CT(-lambdaS + sys.rT*gammaSdot,SS*(mu.*lambdaN));
    error = norm(lambdaN - lambdaN_old) + norm(lambdaS - lambdaS_old);
    converged = error < tol;
    k = k + 1;
end
if ~converged, disp('not converged contact acceleration problem'), end
gammaNdot = (GNN+GNL*QL)*lambdaN + GNS*lambdaS + WN'*(M\h) + nuN;
gammaSdot = (GSN+GSL*QL)*lambdaN + GSS*lambdaS + WS'*(M\h) + nuS;

% index sets of persistent contacts
IP = IN(-lambdaN + sys.rN*gammaNdot<=0); 
IPS = intersect(IP,IS(abs(-lambdaS + sys.rN*gammaSdot)<SS*(mu.*lambdaN)));
IPLplus = intersect(IP,...
    union(ILplus,IS(-lambdaS + sys.rT*gammaSdot>=SS*(mu.*lambdaN))));
IPLmin = intersect(IP,...
    union(ILmin,IS(-lambdaS + sys.rT*gammaSdot<=-SS*(mu.*lambdaN))));

udot = M\(h + (WN + WL*QL)*lambdaN + WS*lambdaS);
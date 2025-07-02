function uplus = impactproblem(sys,t,q,umin,IC,tol,maxiter)
M = sys.M(t,q);
[WN,WT,chiN,chiT] = sys.Wchi(t,q,IC);
[mu,eN,eT] = sys.contactpar(IC);
EN = diag(eN);
ET = diag(eT);
I = eye(length(IC));
GNN = WN'*(M\WN); GNT = WN'*(M\WT); % Delassus submatrices
GTN = WT'*(M\WN); GTT = WT'*(M\WT);
gammaNmin = WN'*umin + chiN;
gammaTmin = WT'*umin + chiT;
converged = 0;
LambdaN = zeros(length(IC),1);
LambdaT = zeros(length(IC),1);
k = 0;
while ~converged && k < maxiter
    LambdaN_old = LambdaN;
    LambdaT_old = LambdaT;
    xiN = GNN*LambdaN + GNT*LambdaT + (EN + I)*gammaNmin;
    xiT = GTN*LambdaN + GTT*LambdaT + (ET + I)*gammaTmin;
    LambdaN = -prox_CN(-LambdaN + sys.rN*xiN);
    LambdaT = -prox_CT(-LambdaT + sys.rT*xiT,mu.*LambdaN);
    error = norm(LambdaN - LambdaN_old) + norm(LambdaT - LambdaT_old);
    converged = error < tol;
    k = k + 1;
end
if ~converged, disp('not converged impact problem'), end
uplus = umin + M\(WN*LambdaN + WT*LambdaT);

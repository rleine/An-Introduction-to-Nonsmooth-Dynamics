function varargout = eventdriven_ode(t,x,flag,sys,IP,IPS,IPLplus,IPLmin)
% eventdriven_ode.m ode_file for eventdriven.m
% Remco Leine, June 2024

q = x(1:sys.ndof,1);
u = x(sys.ndof+1:end,1);
M = sys.M(t,q);
h = sys.h(t,q,u);
n = length(sys.I);

if ~isempty(IP)

    % selection matrices
    S = eye(n);
    SLplus = S(IPLplus,IP);  % forward slipping contacts
    SLmin = S(IPLmin,IP);    % backward slipping contacts
    SS = S(IPS,IP);          % sticking contacts

    [WN,WT,chiN,chiT] = sys.Wchi(t,q,IP); % generalized force directions
    [nuN,nuT] = sys.nu(t,q,u,IP);
    WS = WT*SS';
    WLplus = WT*SLplus';
    WLmin = WT*SLmin';
    WL = [WLplus WLmin];
    nuS = SS*nuT;

    mu = sys.contactpar(IP); % friction coefficients
    QL = [- SLplus*diag(mu); 
            SLmin*diag(mu)]; 

    % Delassus submatrices
    GNN = WN'*(M\WN); GNL = WN'*(M\WL); GNS = WN'*(M\WS);
    GSN = WS'*(M\WN); GSL = WS'*(M\WL); GSS = WS'*(M\WS);
    
    A = [GNN+GNL*QL GNS;
         GSN+GSL*QL GSS];

    b = [WN'*(M\h) + nuN ;
         WS'*(M\h) + nuS ];

    x = -A\b;

    lambdaN = x(1:length(IP),1);
    lambdaS = x(length(IP)+1:end,1);
    lambdaL = QL*lambdaN;

end

switch flag
case ''                          % Return dx/dt = f(t,x).
    varargout{1} = f;
case 'events'                    % Return [value,isterminal,direction]
   [varargout{1:3}] = events;
end


    function dxdt = f
    if ~isempty(IP)  
       udot = M\(h + WN*lambdaN + WS*lambdaS + WL*lambdaL);
    else 
       udot = M\(h);
    end
    dxdt = [u;udot];
    end


    function [value,isterminal,direction] = events
    value = ones(4*n,1);           % default no check
    gN = sys.contactdistance(t,q,sys.I);
    Iopen = setdiff(sys.I,IP);     % index set of open contacts
    value(Iopen,1) = gN(Iopen);    % check for collision at open contacts
    
    if ~isempty(IP)
        value(n+IP,1) = lambdaN;   % check for detachment when closed 
        gammaT = WT'*u + chiT; 

        % check for vanishing sliding velocity during
        value(2*n+IPLplus,1) = SLplus*gammaT; % forward slip 
        value(2*n+IPLmin,1) = -SLmin*gammaT;  % backward slip
        
        % check for saturation of the friction force during stick
        value(3*n+IPS,1) = -abs(lambdaS) + SS*diag(mu)*lambdaN;

    end

    isterminal = ones(4*n,1);
    direction = -ones(4*n,1);
    end

end
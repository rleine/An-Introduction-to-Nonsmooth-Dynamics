function [t,q,u] = eventdriven(sys,tb,te,q0,u0,...
    reltol,abstol,tol,maxiter,trace,tol_index_set)
% eventdriven.m : event driven simulation, Remco Leine, June 2024
if nargin<6, reltol = 1e-8; end
if nargin<7, abstol = 1e-13; end
if nargin<8, tol = reltol/100; end
if nargin<9, maxiter = 1000; end
if nargin<10, trace = 0; end
if nargin<11, tol_index_set = tol;end

q = []; u = []; t = [];
t0 = tb;

% integrator settings
ode_file = 'eventdriven_ode';
opts_ode = odeset('RelTol',reltol,'AbsTol',abstol,'Events','on','Refine',1);

n = length(sys.I); %number of contacts

while t0<te
    
    [WN,WT,chiN,chiT] = sys.Wchi(t0,q0,sys.I); 
    gN = sys.contactdistance(t0,q0,sys.I);
    gammaN = WN'*u0 + chiN;
    gammaT = WT'*u0 + chiT;

    IC = find( gN<=0 );
    IN = intersect(IC,find( gammaN<tol_index_set));
    IS = intersect(IN,find( abs(gammaT)<tol_index_set ));
    ILplus = intersect(IN,find( gammaT>tol_index_set ));
    ILmin = intersect(IN,find( gammaT<-tol_index_set ));
    
    % solve contact acceleration problem to find new index sets
    [IP,IPS,IPLplus,IPLmin] = contactaccelerationproblem(sys,t0,q0,u0,...
        IN,IS,ILplus,ILmin,tol,maxiter);

    % non-impulsive phase  
    [T,X,TE,XE,IE] = ode45(ode_file,[t0 te],[q0;u0],opts_ode,sys,...
        IP,IPS,IPLplus,IPLmin);
   
    if trace, report_event(IE); end
 
    t = [t;T];              % append non-impulsive phase to solution
    q = [q X(:,1:sys.ndof)'];
    u = [u X(:,sys.ndof+1:end)'];
 
    t0 = t(end);
    q0 = q(:,end);
    u0 = u(:,end);

    if ismember(IE,1:n)    % collision => solve impact problem
        gN = sys.contactdistance(t0,q0,sys.I);
        IC = find( gN<=0 );
        u0 = impactproblem(sys,t0,q0,u0,IC,tol,maxiter); 
    end
end

    function event = report_event(IE)
        if isempty(IE)
            event = 'final time reached';
        else
            e_types = {'collision', 'detachment', ...
                'zero sliding velocity', 'friction saturation'};
            contact_num = mod(IE(1)-1,n)+1;
            idx = ceil(IE(1)/n);
            event = sprintf('%s at contact %d',e_types{idx},contact_num);
        end
        fprintf(['t0 = %g, IP = %g, IPS = %g, IPLplus = %g, ' ...
                 'IPLmin = %g, t_end = %g, %s\n'], ...
                 t0, IP, IPS, IPLplus, IPLmin, T(end), event);
    end

end



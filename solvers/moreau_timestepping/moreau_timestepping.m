function [t,q,u] = moreau_timestepping(sys,t0,te,q0,u0,N,tol,maxiter)
% Moreau Timestepping Scheme for systems with 2D friction
% 
% input:  sys        system file in the format of class_sys
%         t0, te     begin and end time
%         q0, u0     initial condition 
%         N          number of discretization points
%         rtol, atol relative/absolute tolerance for fixed point iteration
%         maxiter    maximal number of iterations
%
% output: t          time
%         q          position
%         u          velocity
%
% Remco Leine, May 2025

n = length(sys.I); 
[mu,eN,eT] = deal(sys.mu,sys.eN,sys.eT);

dt = (te-t0)/(N-1);                    %timestep
t = t0:dt:te;                          %time vector

[q,u] = deal(zeros(sys.ndof,N));           %memory allocation 
[PN,PT] = deal(zeros(n,1));  
q(:,1) = q0; u(:,1) = u0;               %initial condition


for i = 1:(N-1)
 
    % start integration step
    qB = q(:,i); 
    uB = u(:,i); 
    tB = t(i);
 
    % calculate midpoint
    qM = qB + dt/2*uB; tM = tB + dt/2;
    hM = sys.h(tM,qM,uB); 
    MM = sys.M(tM,qM);
 
    % calculate index set of closed contacts
    gNM = sys.contactdistance(tM,qM,sys.I);         
    IC = find(gNM<=0)';
 
    % solve contact problem and finish time-step
    if ~isempty(IC)
        [WN,WT,chiN,chiT] = sys.Wchi(tM,qM,IC); 
        gammaNB = WN'*uB + chiN; % contact velocities at tB in IC
        gammaTB = WT'*uB + chiT;
        converged = 0; 
        k = 0;
        while (~converged && k < maxiter) % fixed point iteration
            uE = uB  + MM\(hM*dt + WN*PN(IC) + WT*PT(IC));
            gammaNE = WN'*uE + chiN; % contact velocities at tE in IC
            gammaTE = WT'*uE + chiT;    
            xiN = gammaNE + eN(IC).*gammaNB;
            xiT = gammaTE + eT(IC).*gammaTB;
            PN_new = -prox_CN(-PN(IC) + sys.rN*xiN);
            PT_new = -prox_CT(-PT(IC) + sys.rT*xiT,mu(IC).*PN_new);
            error = max(abs(PN(IC)-PN_new)) + max(abs(PT(IC)-PT_new)); 
            PN(IC) = PN_new;
            PT(IC) = PT_new;
            converged = (error < tol);
            k = k + 1;
        end
        if converged
            uE = uB  + MM\(hM*dt + WN*PN(IC) + WT*PT(IC));
        else    
            fprintf('not converged, t = %g, error = %g\n', tB, error)
        end
    else  %if all constacts are open
        uE = uB  + MM\(hM*dt);
    end
 
    qE = qM + uE*dt/2;

    % update
    q(:,i+1) = qE;
    u(:,i+1) = uE;
end     
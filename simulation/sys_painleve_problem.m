classdef sys_painleve_problem < class_sys

    properties
        m = 1;
        s = 1;
        JS 
        gravity = 10;
        mu
        eN = 0;
        eT = 0;
    end

    methods

        function obj = sys_painleve_problem(mu)
            obj.mu = mu;
            obj.rN = 0.5;
            obj.rT = 0.3;
            obj.ndof = 3;
            obj.I = 1;
            obj.JS = 1/3*obj.m*obj.s^2;
            obj.M_const = diag([obj.m,obj.m,obj.JS]);
            obj.h_const = [ 0;
                           -obj.m*obj.gravity;
                            0];
        end    

        function [WN,WT,chiN,chiT] = Wchi(obj,t,q,I)
            phi = q(3);
     
            WN = [ 0; 1; -obj.s*cos(phi)];
            WT = [ 1; 0; -obj.s*sin(phi)];
     
            WN = WN(:,I);
            WT = WT(:,I);
            chiN = zeros(length(I),1);
            chiT = zeros(length(I),1);
        end


        function [nuN,nuT] = nu(obj,t,q,u,I)
            phi = q(3);
            phidot = u(3);
    
            nuN = obj.s*sin(phi)*phidot^2;
            nuT = -obj.s*cos(phi)*phidot^2;
           
            nuN = nuN(I);
            nuT = nuT(I);
        end

        function gN = contactdistance(obj,t,q,I)
            phi = q(3);
            y = q(2);
            gN = y - obj.s*sin(phi);
            gN = gN(I);
        end    

    end

end


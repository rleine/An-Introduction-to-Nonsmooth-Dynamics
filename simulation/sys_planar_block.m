classdef sys_planar_block < class_sys
% Planar block system 

    properties
        m = 1;
        a
        b
        gravity = 9.81;
        mu
        eN = [0;0];
        eT = [0;0];
    end

    methods

        function obj = sys_planar_block(a,b,k,alpha,mu)
            obj.a = a;
            obj.b = b;
            obj.rN = 1e-2;
            obj.rT = 1e-2;
            obj.ndof = 3;
            obj.I = [1,2];
            obj.mu = [mu;mu];
            J_S = k^2*obj.m;
            obj.M_const = diag([obj.m,obj.m,J_S]);
            obj.h_const = [ obj.m*obj.gravity*sin(alpha);
                           -obj.m*obj.gravity*cos(alpha);
                            0];
        end    

        function [WN,WT,chiN,chiT] = Wchi(obj,t,q,I)
            phi = q(3);
     
            WN(:,1) = [ 0; 1;-obj.a*cos(phi)+obj.b*sin(phi)];
            WN(:,2) = [ 0; 1; obj.a*cos(phi)+obj.b*sin(phi)];
            WT(:,1) = [ 1; 0; obj.a*sin(phi)+obj.b*cos(phi)];
            WT(:,2) = [ 1; 0;-obj.a*sin(phi)+obj.b*cos(phi)];
     
            WN = WN(:,I);
            WT = WT(:,I);
            chiN = zeros(length(I),1);
            chiT = zeros(length(I),1);
        end


        function [nuN,nuT] = nu(obj,t,q,u,I)
            phi = q(3);
            phidot = u(3);
    
            nuN = [(obj.a*sin(phi) + obj.b*cos(phi))*phidot^2;
                   (-obj.a*sin(phi) + obj.b*cos(phi))*phidot^2];

            nuT = [(obj.a*cos(phi) - obj.b*sin(phi))*phidot^2;
                   (-obj.a*cos(phi) - obj.b*sin(phi))*phidot^2];
           
            nuN = nuN(I);
            nuT = nuT(I);
        end

        function gN = contactdistance(obj,t,q,I)
            phi = q(3);
            y = q(2);
            gN1 = y-obj.a*sin(phi)-obj.b*cos(phi);
            gN2 = y+obj.a*sin(phi)-obj.b*cos(phi);
            gN =[gN1;gN2];
            gN = gN(I);
        end    

    end

end


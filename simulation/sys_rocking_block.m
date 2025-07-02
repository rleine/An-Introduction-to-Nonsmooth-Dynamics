classdef sys_rocking_block < class_sys
% Rocking block system 

    properties
        m = 1;
        a
        b
        gravity = 9.81;
        mu = [0.3;0.3];
        eN = [0;0];
        eT = [0;0];
    end

    methods

        function obj = sys_rocking_block(a,b)
            obj.a = a;
            obj.b = b;
            obj.rN = 1e-1;
            obj.rT = 1e-1;
            obj.ndof = 3;
            obj.I = [1,2];
            J = 1/12*obj.m*(4*obj.a^2+4*obj.b^2);
            obj.M_const = diag([obj.m,obj.m,J]);
            obj.h_const = [ 0;
                           -obj.m*obj.gravity;
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


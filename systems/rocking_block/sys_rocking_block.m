classdef sys_rocking_block < class_sys

    properties
        a
        b
    end

    methods

        function obj = sys_rocking_block(a,b,m,J,gravity,mu,r)
            obj.a = a;
            obj.b = b;
            obj.mu = [mu;mu];
            obj.eN = [0;0];
            obj.eT = [0;0];
            obj.rN = r;
            obj.rT = r;
            obj.ndof = 3;
            obj.I = [1,2];
            obj.M_const = diag([m,m,J]);
            obj.h_const = [ 0;
                           -m*gravity;
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
            [chiN,chiT] = deal(zeros(length(I),1));
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


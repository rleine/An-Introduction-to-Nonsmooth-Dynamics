classdef sys_block_on_belt < class_sys

    properties
        m = 1
        k = 1
        gravity = 1
        vdr = 0.2
        delta = 3
    end

    methods

        function obj = sys_block_on_belt
            obj.eN = 0;
            obj.eT = 0;
            obj.mu = 1;
            obj.rN = 0.2;
            obj.rT = 0.2;
            obj.ndof = 2;
            obj.I = 1;
            obj.M_const = diag([obj.m obj.m]);
            obj.WN = [0; 1];  
            obj.WT = [1; 0];
        end    

        function h = h(obj,~,q,u)
            x = q(1);
            xdot = u(1);
            gammaT = xdot - obj.vdr;
            mg = obj.m*obj.gravity;
            FD = obj.mu*mg*obj.delta*gammaT/(1+obj.delta*abs(gammaT));
            h = [FD - obj.k*x; -mg];
        end

        function [WN,WT,chiN,chiT] = Wchi(obj,~,~,I)
            WN = obj.WN(:,I);
            WT = obj.WT(:,I);
            chiN = zeros(length(I),1);
            chiT = -obj.vdr*ones(length(I),1);
        end

        function gN = contactdistance(obj,t,q,I)
            y = q(2);
            gN = y;
            gN = gN(I);
        end    

    end

end
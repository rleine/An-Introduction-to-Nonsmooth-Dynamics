classdef sys_planar_bouncing_ball < class_sys
% planar bouncing ball, Example 3.5

    properties
        R 
    end

    methods

        function obj = sys_planar_bouncing_ball(m,JS,R,gravity,mu,eN,eT)
            obj.rN = 0.2;
            obj.rT = 0.2;
            obj.ndof = 3;
            obj.I = 1;
            obj.R = R;
            obj.mu = mu;
            obj.eN = eN;
            obj.eT = eT;
            obj.M_const = diag([m m JS]);
            obj.h_const = [0; -m*gravity; 0];
            obj.WN = [0; 1; 0];
            obj.WT = [1; 0; obj.R];
        end    

        function gN = contactdistance(obj,t,q,I)
            y = q(2);
            gN = y - obj.R;
            gN = gN(I);
        end    

    end

end
classdef sys_super_ball < class_sys

    properties
        R 
        l
    end

    methods

        function obj = sys_super_ball(m,JS,R,l,gravity,mu,eN,eT)
            obj.rN = 0.2;
            obj.rT = 0.2;
            obj.ndof = 3;
            obj.I = [1 2];
            obj.R = R;
            obj.l = l;
            obj.mu = [mu; mu];
            obj.eN = [eN; eN];
            obj.eT = [eT; eT];
            obj.M_const = diag([m m JS]);
            obj.h_const = [0; -m*gravity; 0];
            obj.WN = [0  0; 
                      1 -1;
                      0  0];
            obj.WT = [1 -1; 
                      0  0; 
                      R  R];
        end    

        function gN = contactdistance(obj,t,q,I)
            y = q(2);
            gN = [y - obj.R;
                  obj.l - y - obj.R];
            gN = gN(I);
        end    

    end

end
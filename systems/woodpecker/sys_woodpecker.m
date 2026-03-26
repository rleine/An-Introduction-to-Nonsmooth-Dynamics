classdef sys_woodpecker < class_sys
% Glocker's model of the Woodpecker Toy 

    properties
        mM = 0.0003
        JM = 5.0e-9
        mS = 0.0045
        JS = 7.0e-7
        lM = 0.010
        lG = 0.015
        lS = 0.0201
        rM = 0.0031
        r0 = 0.0025
        hM = 0.0058
        hS = 0.02
        gravity = 9.81
        cphi = 0.0056
    end

    methods

        function obj = sys_woodpecker
            obj.mu = [0.3;0.3;0.3];
            obj.eN = [0.5;0;0];
            obj.eT = [0;0;0];
            obj.rN = 1e-4;
            obj.rT = 1e-4;
            obj.ndof = 3;
            obj.I = [1,2,3];
            obj.WN(:,1) = [ 0; 0;-obj.hS];
            obj.WN(:,2) = [ 0; obj.hM; 0];
            obj.WN(:,3) = [ 0;-obj.hM; 0];
            obj.WT(:,1) = [ 1; obj.lM; obj.lG-obj.lS];
            obj.WT(:,2) = [ 1; obj.rM; 0];
            obj.WT(:,3) = [ 1; obj.rM; 0];
            mS = obj.mS; mM = obj.mM;
            lM = obj.lM; lG = obj.lG;
            JM = obj.JM; JS = obj.JS;
            obj.M_const = [ mS+mM  mS*lM      mS*lG       ;
                            mS*lM  JM+mS*lM^2 mS*lM*lG    ;
                            mS*lG  mS*lM*lG   JS+mS*lG^2 ];
        end    


        function h = h(obj,t,q,u)
            phiM = q(2);
            phiS = q(3);
            h = [ -(obj.mS+obj.mM)*obj.gravity;
                 -obj.cphi*(phiM-phiS)-obj.mS*obj.gravity*obj.lM;
                 -obj.cphi*(phiS-phiM)-obj.mS*obj.gravity*obj.lG];
        end

        function gN = contactdistance(obj,t,q,I)
            phiM = q(2);
            phiS = q(3);
            gN1 = (obj.lM+obj.lG-obj.lS-obj.r0) - obj.hS*phiS;
            gN2 = (obj.rM-obj.r0) + obj.hM*phiM;
            gN3 = (obj.rM-obj.r0) - obj.hM*phiM;
            gN = [gN1;gN2;gN3];
            gN = gN(I);
        end    

    end

end
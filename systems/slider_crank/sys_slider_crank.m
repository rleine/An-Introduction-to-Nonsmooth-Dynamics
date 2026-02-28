classdef sys_slider_crank < class_sys

    properties     
        l1 = 0.153 
        l2 = 0.306
        a = 0.05
        b = 0.025
        c = 0.001
        d = 0.052 
        m1 = 0.038
        m2 = 0.038
        m3 = 0.076
        J1 = 7.4e-5
        J2 = 5.9e-4
        J3 = 2.7e-6
        gravity = 9.81;
    end

    methods

        function obj = sys_slider_crank(mu,eN,eT)
            obj.rN = 0.001;
            obj.rT = 0.001;
            obj.ndof = 3;
            obj.I = [1 2 3 4];
            obj.mu = [mu; mu; mu; mu];
            obj.eN = [eN; eN; eN; eN];
            obj.eT = [eT; eT; eT; eT];
        end    

        function M = M(obj,~,q)
            m1 = obj.m1; m2 = obj.m2; m3 = obj.m3;
            J1 = obj.J1; J2 = obj.J2; J3 = obj.J3;
            l1 = obj.l1; l2 = obj.l2;
            M = zeros(3,3);
            M(1,1) = J1 + (m1/4 + m2 + m3)*l1^2;
            M(1,2) = (m2/2 + m3)*l1*l2*cos(q(2) - q(1));
            M(2,1) = M(1,2);
            M(2,2) = J2 + (m2/4 + m3)*l2^2;
            M(3,3) = J3;
        end

        function h = h(obj,~,q,u)
            m1 = obj.m1; m2 = obj.m2; m3 = obj.m3;
            l1 = obj.l1; l2 = obj.l2;
            g = obj.gravity;
            h_gyro = [ (m2/2 + m3)*l1*l2*sin(q(2)-q(1))*u(2)^2;
                      -(m2/2 + m3)*l1*l2*sin(q(2)-q(1))*u(1)^2;
                       0];
            h_grav = [ - (m1/2 + m2 + m3)*g*l1*cos(q(1));
                       - (m2/2 + m3)*g*l2*cos(q(2));
                       0];
            h = h_gyro + h_grav;
        end

        function [WN,WT,chiN,chiT] = Wchi(obj,~,q,I)
            l1 = obj.l1; l2 = obj.l2;
            a = obj.a;    b = obj.b;
            sn = sin(q); cs = cos(q);
            wN1 = [-l1*cs(1); -l2*cs(2);  a*cs(3) + b*sn(3)];
            wN2 = [-l1*cs(1); -l2*cs(2); -a*cs(3) + b*sn(3)];
            wN3 = [ l1*cs(1);  l2*cs(2); -a*cs(3) + b*sn(3)];
            wN4 = [ l1*cs(1);  l2*cs(2);  a*cs(3) + b*sn(3)];   
            wT1 = [-l1*sn(1); -l2*sn(2);  a*sn(3) - b*cs(3)];
            wT2 = [-l1*sn(1); -l2*sn(2); -a*sn(3) - b*cs(3)];
            wT3 = [-l1*sn(1); -l2*sn(2);  a*sn(3) + b*cs(3)];
            wT4 = [-l1*sn(1); -l2*sn(2); -a*sn(3) + b*cs(3)];
            WN = [wN1 wN2 wN3 wN4];
            WT = [wT1 wT2 wT3 wT4];
            WN = WN(:,I);
            WT = WT(:,I);
            [chiN,chiT] = deal(zeros(length(I),1));
        end

        function [nuN,nuT] = nu(obj,~,q,u,I) 
            l1 = obj.l1; l2 = obj.l2;
            a = obj.a;    b = obj.b;
            sn = sin(q); cs = cos(q);
            nuN1 = [l1*sn(1)  l2*sn(2) (-a*sn(3) + b*cs(3))]*u.^2;
            nuN2 = [l1*sn(1)  l2*sn(2) ( a*sn(3) + b*cs(3))]*u.^2;
            nuN3 = [l1*sn(1) -l2*sn(2) ( a*sn(3) + b*cs(3))]*u.^2;
            nuN4 = [l1*sn(1) -l2*sn(2) (-a*sn(3) + b*cs(3))]*u.^2;
            nuN = [nuN1;nuN2;nuN3;nuN4];
            nuT1 = [-l1*cs(1) -l2*cs(2) ( a*cs(3) + b*sn(3))]*u.^2;
            nuT2 = [-l1*cs(1) -l2*cs(2) (-a*cs(3) + b*sn(3))]*u.^2;
            nuT3 = [-l1*cs(1) -l2*cs(2) ( a*cs(3) - b*sn(3))]*u.^2;
            nuT4 = [-l1*cs(1) -l2*cs(2) (-a*cs(3) - b*sn(3))]*u.^2;
            nuT = [nuT1;nuT2;nuT3;nuT4];        
            nuN = nuN(I);
            nuT = nuT(I);
        end

        function gN = contactdistance(obj,~,q,I)
            sn = sin(q); cs = cos(q);
            l1 = obj.l1; l2 = obj.l2;
            a = obj.a;   b = obj.b;  d = obj.d;
            gN = [d/2 - l1*sn(1) - l2*sn(2) + a*sn(3) - b*cs(3);
                  d/2 - l1*sn(1) - l2*sn(2) - a*sn(3) - b*cs(3);
                  d/2 + l1*sn(1) + l2*sn(2) - a*sn(3) - b*cs(3);
                  d/2 + l1*sn(1) + l2*sn(2) + a*sn(3) - b*cs(3)];
            gN = gN(I);
        end    

    end

end
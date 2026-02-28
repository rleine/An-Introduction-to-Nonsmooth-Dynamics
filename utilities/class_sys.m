classdef class_sys
% Base class of a mechanical system with frictional unilateral constraints

    properties
        M_const
        h_const
        WN
        WT
        I
        ndof
        mu 
        eN
        eT
        rN
        rT
    end

    methods
       
        function M = M(obj,~,~)
            M = obj.M_const;
        end

        function h = h(obj,~,~,~)
            h = obj.h_const;
        end

        function [WN,WT,chiN,chiT] = Wchi(obj,~,~,I)
            WN = obj.WN(:,I);
            WT = obj.WT(:,I);
            [chiN,chiT] = deal(zeros(length(I),1));
        end

        function [nuN,nuT] = nu(obj,~,~,~,I)
            [nuN,nuT] = deal(zeros(length(I),1));
        end

        function [mu,eN,eT] = contactpar(obj,I)
            mu = obj.mu(I);
            eN = obj.eN(I);
            eT = obj.eT(I);
        end    

        function gN = constactdistance(obj,t,q,I)
 
        end    

        function gNall = getgNall(obj,t,q) %for post-processing
            gNall = zeros(length(obj.I),length(t));
            for i = 1:length(t)
                 gNall(:,i) = obj.contactdistance(t(i),q(:,i),obj.I);
            end     
        end    
    end

end
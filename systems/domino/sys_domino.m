classdef sys_domino < class_sys

    properties
        n
        a
        b
    end

    methods

        function obj = sys_domino(n,a,b,m,J,gravity,muG,muD,eG,eD,r)
            obj.n = n;
            obj.a = a;
            obj.b = b;
            obj.mu = [muG*ones(3*n,1);muD*ones(2*n-2,1)];
            obj.eN = [eG*ones(3*n,1);eD*ones(2*n-2,1)];
            obj.eT = zeros(5*n-2,1);
            obj.rN = r;
            obj.rT = r;
            obj.ndof = 3*n;
            obj.I = 1:(5*n-2);
            obj.M_const = kron(eye(n),diag([m,m,J]));
            obj.h_const = kron(ones(n,1),[ 0; -m*gravity; 0]);
        end    

        function [WN,WT,chiN,chiT] = Wchi(obj,~,q,I)
            [W_NL,W_NR,W_NT,W_TL,W_TR,W_TT] = deal(zeros(obj.ndof,obj.n));
            [W_NB,W_ND,W_TB,W_TD] = deal(zeros(obj.ndof,obj.n-1));

            x = q(1:3:end);
            y = q(2:3:end);
            phi = q(3:3:end);
            c = cos(phi); s = sin(phi);

            xB = x - obj.a*c + obj.b*s;
            yB = y - obj.a*s - obj.b*c;

            xD = x + obj.a*c - obj.b*s;
            yD = y + obj.a*s + obj.b*c;

            for i=1:obj.n
                C_i = zeros(3, obj.ndof); C_i(:, 3*(i-1)+1 : 3*i) = eye(3);
             
                dgNLidqi = [ 0 1 -obj.a*c(i)+obj.b*s(i)];
                dgNRidqi = [ 0 1  obj.a*c(i)+obj.b*s(i)];
                dgNTidqi = [ 0 1  obj.a*c(i)-obj.b*s(i)];
                dgTLidqi = [ 1 0  obj.a*s(i)+obj.b*c(i)];
                dgTRidqi = [ 1 0 -obj.a*s(i)+obj.b*c(i)];
                dgTTidqi = [ 1 0 -obj.a*s(i)-obj.b*c(i)];
                
                W_NL(:,i) = (dgNLidqi*C_i)';
                W_NR(:,i) = (dgNRidqi*C_i)';
                W_NT(:,i) = (dgNTidqi*C_i)';
                W_TL(:,i) = (dgTLidqi*C_i)';
                W_TR(:,i) = (dgTRidqi*C_i)';
                W_TT(:,i) = (dgTTidqi*C_i)';
                
                if i<obj.n && any(ismember([3*obj.n+i, 4*obj.n-1+i], I))

                    j = i+1;
                    C_j = zeros(3, obj.ndof); C_j(:, 3*(j-1)+1 : 3*j) = eye(3);

                    dgNBidqi = [ -c(i) -s(i) -s(i)*(xB(j)-x(i))+c(i)*(yB(j)-y(i))];
                    dgNBidqj = [  c(i)  s(i)  s(i)*(xB(j)-x(j))-c(i)*(yB(j)-y(j))];
                    dgTBidqi = [  s(i) -c(i) -c(i)*(xB(j)-x(i))-s(i)*(yB(j)-y(i))];
                    dgTBidqj = [ -s(i)  c(i)  c(i)*(xB(j)-x(j))+s(i)*(yB(j)-y(j))];

                    dgNDidqi = [ -c(j) -s(j) -s(j)*(xD(i)-x(i))+c(j)*(yD(i)-y(i))];
                    dgNDidqj = [  c(j)  s(j)  s(j)*(xD(i)-x(j))-c(j)*(yD(i)-y(j))];
                    dgTDidqi = [  s(j) -c(j) -c(j)*(xD(i)-x(i))-s(j)*(yD(i)-y(i))];
                    dgTDidqj = [ -s(j)  c(j)  c(j)*(xD(i)-x(j))+s(j)*(yD(i)-y(j))];

                    W_NB(:,i) = (dgNBidqi*C_i)' + (dgNBidqj*C_j)';
                    W_TB(:,i) = (dgTBidqi*C_i)' + (dgTBidqj*C_j)';
                    W_ND(:,i) = (dgNDidqi*C_i)' + (dgNDidqj*C_j)';
                    W_TD(:,i) = (dgTDidqi*C_i)' + (dgTDidqj*C_j)';
                end
            end
            
            WN = [W_NL W_NR W_NT W_NB W_ND];
            WT = [W_TL W_TR W_TT W_TB W_TD];
     
            WN = WN(:,I);
            WT = WT(:,I);
            chiN = zeros(length(I),1);
            chiT = zeros(length(I),1);

        end

        function gN = contactdistance(obj,~,q,I)

            x = q(1:3:end);
            y = q(2:3:end);
            phi = q(3:3:end);
            c = cos(phi); s = sin(phi);

            xB = x - obj.a*c + obj.b*s;
            yB = y - obj.a*s - obj.b*c;

            xD = x + obj.a*c - obj.b*s;
            yD = y + obj.a*s + obj.b*c;

            g_NL = yB;
            g_NR = y + obj.a*s - obj.b*c;
            g_NT = yD;

            ind_i = 1:obj.n-1;
            ind_j = 2:obj.n;
            g_NB =  c(ind_i).*(xB(ind_j) - x(ind_i)) + s(ind_i).*(yB(ind_j) - y(ind_i)) - obj.a;
            g_ND = -c(ind_j).*(xD(ind_i) - x(ind_j)) - s(ind_j).*(yD(ind_i) - y(ind_j)) - obj.a;

            gN = [g_NL;g_NR;g_NT;g_NB;g_ND];
            gN = gN(I);
        end    

    end

end

% Description: Lie is a static class that conatins basic
%              routines and properties to manuplate SE(n)
%              and its direct products.

% 2017  Girum G. Demisse, girumdemisse@gmail.com/girum.demisse@uni.lu
%       Computer vision team, University of Luxembourg.
%------------------------------------------------

classdef Lie
    
    methods(Static)
         
        
        function arr =  arrange(rep,val)
            % Description: Performs permutation of a direct product
            %              of rigid-transformation matrices.
            % INPUT      : rep  - (nxnxm) is m-direct product of
            %                     n-dimensional rigid transformation matrices.
            %            : val  - m-dimensional vector representation
            %                     of a permutation.
            % OUTPUT     : arr  - (nxnxm) permuted "rep" according to "val"
            %                     rigid-transformations.
            
            [~,c,r] = size(rep);
            arr = zeros(c,c,r);
            for k=1:r
                T = rep(:,:,val(k));
                arr(:,:,k) = T;
            end
        end
        
        
        
        function [res2, vec] = pull_pro_lie(G1,G2)
            % Description: Pulls the second argument to the tangent space of the first.
            %
            % INPUT      : G1  - (nxnxm) n-dimensional rigid transformation matrix.
            %                    Note: G1 can be empty. If  so G2 will be pulled
            %                    to the tangent space of the idenitity(Lie algebra).
            %
            %            : G2  - (nxnxm) n-dimensional rigid transformation matrix.
            % OUTPUT     : res - (nxnxm) G2 pulled to tangent space of G1.
            %            : vec - (n^2x1) res in vector form.
            
            [n,~,m] = size(G2);
            res2 = [];
            vec = [];
            %dim = n^2;
            if isempty(G1)
                for i=1:m
                    [temp,VecR] = Lie.pull_lie(eye(n),G2(:,:,i));
                    res2 = cat(3,res2, temp);
                    vec = [vec;VecR];
                end
            else
                for i=1:m
                    [temp,VecR] = Lie.pull_lie(G1(:,:,i),G2(:,:,i));
                    res2 = cat(3,res2, temp);
                    vec = [vec;VecR];
                end
            end
        end
        
        
        
        function [res,vec] = pull_lie(MAT1,MAT2)
            % Description: Pulls the second argument to the tangent space of the first.
            %
            % INPUT      : MAT1  - (nxn) n-dimensional rigid transformation matrix.
            %            : MAT2  - (nxn) n-dimensional rigid transformation matrix.
            
            % OUTPUT     : res   - (nxn) MAT2 pulled to tangent space of MAT1
            %            : vec   - (n^2x1) res in vector form.
            
            [n,~] = size(MAT1);
            %dim = n^2;
            rot  = real(logm(MAT1(1:n-1,1:n-1)'*MAT2(1:n-1,1:n-1)));
            tra  = MAT2(1:n-1,1:n-1)'*(MAT2(1:n-1,n) - MAT1(1:n-1,n));
            
            res  =[rot,tra; zeros(1,n)];
            %vec = reshape(res,[dim, 1]);
            vec = [rot,tra];
            vec = vec(:);
        end
        
        
        
        function res = group_pro(MAT1,MAT2)
            % Description: Element-wise group product
            %
            % INPUT      : MAT1  - (nxnxm) m n-dimensional rigid transformation
            %                       matrices.
            % OUTPUT     : res   - (nxnxm) element-wise product.
            
            [~,~,m] = size(MAT1);
            
            res = [];
            for i=1:m
                temp = MAT1(:,:,i)*MAT2(:,:,i);
                res = cat(3,res,temp);
            end
        end
        
        
        
        function res = group_invert(MAT1)
            % Description: inverts a direct product of rigid transformations
            %
            % INPUT      : MAT1 - (nxnxm) m n-dimensional rigid transformation
            %                      matrices.
            % OUTPUT     : res  - (nxnxm) element-wise inverse.
            
            m = size(MAT1,3);
            res = [];
            for i=1:m
                temp = inv(MAT1(:,:,i));
                res = cat(3,res,temp);
            end
        end
        
        
        
        function icur = inver(G1)
            % Description: Pair-wise inversion and order flipping of a direct product
            %              of rigid-transformation matrices.
            % INPUT      : G1   - (nxnxm) is m direct product of
            %                      n-dimensional rigid transformation matrices.
            % OUTPUT     : icur - (nxnxm) pair-wise inverse of the input in reverse order.
            
            
            [~,c,t] = size(G1);
            icur = zeros(c,c,t);
            for k=1:t
                l = t-k+1;
                temp = G1(:,:,l);
                icur(:,:,k) = inv(temp);
            end
        end
        
        
        
        function val = geo_dist_SE(MAT1,MAT2)
            % Description: Computes  the geodeisc distance between
            %              two rigid-transformation matrices.
            % INPUT      : MAT1 - (nxn) n-dimensional rigid transformation
            %                      matrice.
            %            : MAT2 - (nxn) n-dimensional rigid transformation
            %                     matrice.
            % OUTPUT     : dist - a scalar distance value.
            
            [~,n] = size(MAT1);
            % weigh each component if necessary
            dM = 1;
            rM = 1;
            tr = n-1;
            rot = rM*logm(MAT1(1:tr,1:tr)'*MAT2(1:tr,1:tr));
            trans  = dM*(MAT1(1:tr,n)-MAT2(1:tr,n));
            val = ((norm(rot))^2 + (norm(trans))^2)^0.5;
        end
        
        
        
        function dist = dist_pro_lie(MAT1,MAT2)
            % Description: Computes  the geodeisc distance between
            %              direct product of rigid-transformation matrices.
            % INPUT      : MAT1   - (nxnxm) m direct product of n-dimensional rigid
            %                       transformation matrices.
            %            : MAT2   - (nxnxm) m direct product of n-dimensional rigid
            %                       transformation matrices.
            %
            % OUTPUT     : dist   - a scalar distance value.
            
            [~,~,m] = size(MAT1);
            dist =  0;
            for j=1:m
                Temp = (Lie.geo_dist_SE(MAT1(:,:,j),MAT2(:,:,j)))^2;
                dist = dist + Temp;
            end
            dist = (dist)^0.5;
        end
        
        
        
        function mean = meanSE(Mat,weight)
            % Description: Computes karcher mean of rigid-transformation
            %              matrices
            % INPUT      : Mat    - (nxnxm) m, n-dimensional rigid
            %                       transformation matrices.
            %            : weight - (mx1) weighting vector. CAUTION sum(weight) should be 1
            % OUTPUT     : mean   -  nxn rigid-transformation matrix.
            
            [r,~,m] = size(Mat);
            rot = Mat(1:r-1,1:r-1,:);
            
            
            
            %--- inital guess of the mean rotation is computed based on
            % Ho jeffery, et.al. "Recursive Karcher Expectation Estimators And Geometric Law
            % of Large Numbers" AISTATS, 2013.
            %-----
            
            temp = rot(:,:,1);
            meanT = weight(1,1)*Mat(1:r-1,r,1);
            
            for j=2:m
                temp = temp*(temp'*rot(:,:,j))^(1/(j+1));
                % arithmetic mean of translation
                meanT = meanT + weight(j,1)*Mat(1:r-1,r,j);
            end
            gmean = real(temp);
            
            
            % free parameters
            
            thresh = 0.00002;        % gradient tolerance: increase for high accuracy
            grad = zeros(r-1,r-1);   % inital gradient value of the karcher equation
            step_size = 0.0005;      % learning rate. CAUTION: high value might lead to jumping back and forth over the optimal value.
            gsearch = 1;             % flag turns off when optimal value is found
            count = 0;               % itteration counter
            stoptime = 100;          % Allowed itteration, if frobenius norm of gradient is NOT smaller than "thresh" by stoptime search HALTS
                                     % and current value will be returned.
            
            % inital guess failure precaution
            if norm(gmean,2) == 0
                gmean = eye(r-1);
            end
            
            % start searching
            while gsearch
                
                for j=1:m
                    grad = grad + weight(j,1)*logm(rot(:,:,j)'*gmean);
                end
                
                val = norm(grad,2);
                
                if val > thresh
                    gmean= gmean * expm(step_size*-1*grad);
                    grad = zeros(r-1,r-1);
                else
                    gsearch = 0;
                end
                % counting the iteration
                if count > stoptime
                    break;
                end
                count = count + 1;
                
            end
            
            % translation mean
            MM = [real(gmean),meanT];
            MM = [MM; zeros(1,r-1),1];
            mean = MM;
        end
    end
end


% Description: Tools is a static class that conatins basic
%              routines and helper functions for Curve class and more.

% 2017  Girum G. Demisse, girumdemisse@gmail.com/girum.demisse@uni.lu
%       Computer vision team, University of Luxembourg.
%------------------------------------------------


classdef Tools
    
    methods(Static)
        
        function [rot,trans] = Trans_Mat(point1,point2)
            % Description: Computes high dimensional rigid-transformation
            %              from point1 to point 2 that preserves world coordinate frame.
            % INPUT      : point1, and point2 - two row vectors (1xn)
            % OUTPUT     : rot                - (nxn) rotation matrix
            %              trans              - (nx1) translation vector
            
            c = size(point1,2);
            CC = [point1;point2];
            CC = CC - repmat(mean(CC),2,1);
            
            % rotation plane estimation
            [~,~,V]=svd(CC);
            
            p1 = (V'*point1')'; p2 = (V'*point2')';
            
            % solve optimal rotation on 2D projection of points
            R = Tools.Opt_rot(p1(1,1:2),p2(1,1:2));
            
            E = eye(c);
            
            E(1:2,1:2) = R;
            rot = E;
            
            % representing rotation matrix with respect to original frame
            rot = V*rot*V';
            
            trans  =  point2' - rot*point1';
        end
        
        
        
        function rot = Opt_rot(p1, p2)
            % Description: computes optimal rotation from p1 to p2
            % INPUT      : p1 and p2 - two row vectors (1xn)
            % OUTPUT     : rot       - (nxn) rotation matrix
            
            c = size(p1,2);
            corr = p1'*p2;
            [U,~,V]=svd(corr);
            E = eye(c);
            E(end,end) = det(V*U');
            rot = V*E*U';
        end
        
        
        
        function [area,comps] = Curve_area(points,varargin)
            % Description: computes area of a PLANAR curve. if the curve is closed 
            %              it uses Green's theorem.
            %
            % INPUT      : points     - ordered set of 2-dimensional points (mx2)
            %              [optional] - set 0 for open curves
            %                         - set 1 for closed curves (default).
            % OUTPUT     : area       - final result (scalar)
            %              comps      - piece-wise area under the curve (mx1)
            
            pointsT = [points;points(1,:)];
            
            if nargin > 1
                if varargin{1} == 0
                  % CAUTION !! integration along the first dimension
                  pointsT = points;
                end
            end
                 
            
            avg_y = 0.5 * (pointsT(1:end-1,2) + pointsT(2:end,2));
            integral = diff(pointsT(:,1)).* avg_y;
            comps = integral;
            
            % direction of point ordering does not matter due to
            % the absolute value.
            area = abs(sum(comps));
            
        end
        
        
        function [integrated,comps] = Curve_length(points,varargin)
            % Description: computes length of a curve
            % INPUT      : points            - ordered set of n-dimensional points (mxn)
            %            : [optional] scalar - 0 to identify open curves
            %                                - different from 0 to
            %                                  identify closed curves (default).
            % OUTPUT     : integrated        - final curve length (scalar)
            %              comps             - piece-wise length (mx1)
            
            if nargin > 1
                % curve length for open curves
                if varargin{1} == 0
                    tangent = diff(points);
                end
            else
                % curve length for closed curves
                tangent = (points - [points(end,:);points(1:end-1,:)]);
            end
            
            comps = sqrt(sum(tangent.^2,2));
            integrated = sum(comps);
        end
        
        
        
        function [m,ind] =  Min_M(MAT)
            % Description: Computes minimum value from an array of reals
            % INPUT      : MAT -(mxn) array of numbers
            % OUTPUT     : m   - the minimum value
            %              ind - index of the solution, [row_index,column_index]
            
            ind = zeros(1,2);
            [l,I] = min(MAT);
            [m,ind(1,2)] = min(l);
            ind(1,1) = I(ind(1,2));
        end
        
        
        
        function points = process(points,varargin)
            % Description: The function samples and normaliz points from it
            %              argument "points".
            % INPUT      : points     - (mxn)set of m input points
            %            : [optional] - scalar to state the sample size.
            % OUTPUT     : points     - sampled and normalized points
            
            
            if nargin > 1
                % number of sample size is specified
                points = curvspace(points,varargin{1});
            else
                % Default: takes input size as sample size
                points = curvspace(points,length(points));
            end
            points = Tools.normal_points(points);
        end
        
        
        function Normal_points = normal_points(points)
            % Description: The function centers and normalizes input data.
            % INPUT      : points         - (mxn) orderd set of m input points in
            %                                n-dimensional space.
            % OUTPUT     : Normals_points - (mxn) centered and normalized points
            
            
            r = size(points,1);
            cen_p = points - repmat(mean(points),r,1);
            Normal_points   = cen_p./sqrt(sum(sum(cen_p.^2,2)));
        end
        

        function [c1_r,c2_r] = DP_sampling(c1,c2,alpha,beta,window_size,type,sample_points)
            % Description: Optimally sample points for correspondance
            %              estimation.
            % INPUT      : c1             - instance of a Curve class.
            %            : c2             - instance of a Curve class, NOTE: |c2.points| == |c1.points|
            %            : alpha          - a scalar constant to weight deformation term.
            %            : beta           - a scalar constant to weight constraint term.
            %            : window_size    - a positive scalar determing the window size of the
            %                               charts.
            %            : type           - = 1 for fixing c1 to uniform and optimally sampling c2.
            %                             - = 0 for sampling both curves optimally.
            %            : sample_points  - scalar value, determins the number of points to be sampled
            %                               sample_points << |c2.points|.
            % OUTPUT     : c1_r           - incase of "type" = 1, uniformly sampled c1.
            %                               incase of "type" = 0, optimally sampled c1.
            %            : c2_r           - optimally sampled  c2.
            
            num_points = size(c1.points,1);
            
            if  sample_points > num_points
                error('The number of sample points exceeds the total number of points !!');
            elseif sample_points <= 0
                error('The number of sample points cannot be zero or negartive !!');
            elseif mod(num_points,sample_points) ~= 0
                warning('The number of sample points does not divide the total number of points. Choose approprately for exact solution.');
            end
            
            
            if type == 1
                
                l = ceil(num_points/sample_points);
                % uniformly sampled c1
                c1_r = c1;
                c1_r.points = Tools.normal_points(c1_r.points(1:l:num_points,:));
                %c1_r.points =  Tools.process(c1.points,sample_points);
                
                c1_r = curve_rep(c1_r);
                
                opt_index = dp_optimal_point_sampling_Single(c1_r,c2,alpha,beta,window_size);
                
                % second curve optimal sampling solution
                c2_r = c2;
                c2_r.points = Tools.normal_points(c2.points(opt_index,:));
                
                c2_r = curve_rep(c2_r);
                
            elseif type == 0
                % sampling both shapes
                opt_index = dp_optimal_point_sampling_Both(c1,c2,alpha,beta,window_size,sample_points);
                
                % first curve optimal sampling solution
                c1_r = c1;
                c1_r.points =  Tools.normal_points(c1.points(opt_index(:,1),:));
                c1_r = curve_rep(c1_r);
                
                % second curve optimal sampling solution
                c2_r = c2;
                c2_r.points =  Tools.normal_points(c2.points(opt_index(:,2),:));
                c2_r = curve_rep(c2_r);
            end
        end
        
        
        
        function plot_opt_paths(cost,index1,index2,index3)
            % Description: Visualizes the optimal solution and the
            %              search space defined by the window size.
            %
            % INPUT      : cost   - (m,n) or (m,m,n) dimensional array of
            %            : index1 - (n,1) vector, solution based on
            %                        deformation only.
            %            : index2 - (n,1) vector, solution based on
            %                        area only.
            %            : index3 - (n,1) vector, solution based on
            %                        constrained objective functional.
            
            
            R = cost;
            sample_size = size(R,3);
            
            figure;
            if sample_size == 1
                % heat plot of cost and solution paths
                imagesc(R); hold on; 
                set(gca,'XAxisLocation','top');
                xllhand = get(gca,'xlabel');  yllhand = get(gca,'ylabel');
                set(xllhand,'string','Index of sampled points: [1,z]','fontsize',15);
                set(yllhand,'string','Index of shape points','fontsize',15);
                colormap jet;
                hold on;
                %solution paths
                plot(index1,'Linewidth',2,'Color','green'); hold on; % deformation
                plot(index2,'Linewidth',2,'color','red');hold on; % area
                plot(index3,'Linewidth',2,'color','yellow'); % constrained
                legend('Deformation ONLY','Area preservation ONLY','Combined');
            else
                % a sliced visualization of the cost space and the solution path
                slices = 6;
                l = ceil(sample_size/slices);
                Sz = 1:l:sample_size;
                h = slice(R,[],[],Sz);
                set(h,'FaceColor','interp','EdgeColor','none'); alpha 0.25;
                contourslice(R,[],[],Sz);
                colormap jet;
                hold on;
                f_inds = (1:1:sample_size)'; % index of the subfunctions
                line(index1(:,1),index1(:,2),f_inds,'DisplayName','Deformation ONLY','Linewidth',2,'Color','green'); hold on; % unconstrained
                line(index2(:,1),index2(:,2),f_inds,'Linewidth',2,'color','red'); hold on; % uniform sampling
                line(index3(:,1),index3(:,2),f_inds,'Linewidth',2,'color','yellow'); hold on; %  constrained
                %legend('','Area preservation ONLY','Combined');
                
                ax = gca;
                ax.BoxStyle = 'full';
                box on; hold on; grid off; hold on;
                set(gca,'XAxisLocation','top');
                xllhand = get(gca,'xlabel');  yllhand = get(gca,'ylabel'); zllhand = get(gca,'zlabel');
                set(xllhand,'string','Curve-1','fontsize',10);
                set(yllhand,'string','Curve-2','fontsize',10);
                set(zllhand,'string','Sub-functions  index [1-z]','fontsize',10); hold on;
                view([-37,32]);
            end
            title('Solution paths');
            hold off;
        end
        
    end
end
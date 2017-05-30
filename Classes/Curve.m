

% Description: Curve is a value class that conatins basic
%              routines and properties to use and manuplate
%              the curved shape representation proposed in [1,2].

% 1) Demisse, G.G, Aouada,D, Ottersten, B. "Similarity Metric For Curved
%    Shapes in Eculidean Space.", IEEE CVPR 2016.
% 2) Demisse, G.G, Aouada,D, Ottersten, B. "Deformation based Curved Shape
%    Representation". IEEE TPAMI 2017.
%
% 2017  Girum G. Demisse, girumdemisse@gmail.com/girum.demisse@uni.lu
%       Computer vision team, University of Luxembourg.
%------------------------------------------------





classdef Curve
    
    properties
        points        % uniformly sampled points
        start_point   % fixed starting point
        type          % type of curve ( 0 - for open, 1- for closed)
        reps          % curve representation
    end
    
    properties(Constant)
        %-- NOTE: A simple 2-dimensional finite (order 4) Abelian group of reflection is defined here.
        %--       We use this definitaion to achive reflection invariance in 2D
        %--       space.
        
        AbelianR = cat(3, [1 0;...
                            0 1], ...
                            [-1 0;...
                            0 1],...
                            [1 0;...
                            0 -1],...
                            [-1 0;...
                            0 -1]);
    end
    
    methods
        
        function ob1 = Curve(data,varargin)
            % Description: Class constructor
            % INPUT      : data       - (mxn) ordered array of m-points in
            %                            n-dimensional space.
            %            : [optional] - type, 0 (open curve) or 1 (closed curve).
            %            : [optional] - start_point, (1xn) fixed starting point
            %                           e.g. of initalization: Curve(points) OR
            %                           Curve(points, type) OR Curve(points,type,start_point)
            % OUTPUT     : o1b        - instance of a Curve class
            
            ob1.points =  data;
            try
                switch nargin
                    case 2
                        if isscalar(varargin{1})
                            ob1.type = varargin{1};
                            l = zeros(1,size(data,2));
                            l(1,1) = 1;
                            ob1.start_point = l;
                        else
                            M = MException('type mismatch','');
                            throw(M);
                        end
                    case 3
                        if isscalar(varargin{1}) && ~isscalar(varargin{2})
                            ob1.type = varargin{1};
                            ob1.start_point = varargin{2};
                        else
                            M = MException('type mismatch','');
                            throw(M);
                        end
                    otherwise
                        % dafault starting point and curve type
                        l = zeros(1,size(data,2));
                        l(1,1) = 1;
                        ob1.start_point = l;
                        ob1.type = 1;
                end
            catch
                error('Make sure you passed the correct parameter values in the right order: data, [type], [start_point]')
            end
            % curve representation
            if ~isempty(ob1.points)
                ob1 = curve_rep(ob1);
            end
        end
        
        
        
        function  ob1 = curve_rep(ob1)
            % Description: computes the representation of a given curve
            % INPUT      : ob1 - instance of curve object
            % OUTPUT     : ob1 - a copy of the input object with initalized "reps" property
            
            if isempty(ob1.points)
                error('Initalize the points property first.');
            end
            
            % fixing the starting point to the refernce
            dist = ob1.points(1,:) - ob1.start_point;
            ob1.points = ob1.points - repmat(dist,length(ob1.points),1);
            
            [r,c] = size(ob1.points);
            
            % closed curve representation
            if ob1.type == 1
                rep = zeros(c+1,c+1,r);
                for i=1:r
                    if i == r
                        [rot,trans] = Tools.Trans_Mat(ob1.points(i,:),ob1.points(1,:));
                        rep(:,:,i)  = [rot,trans;...
                            zeros(1,c),1];
                    else
                        [rot,trans] = Tools.Trans_Mat(ob1.points(i,:),ob1.points(i+1,:));
                        rep(:,:,i)  = [rot,trans;...
                            zeros(1,c),1];
                    end
                end
                % open curve representation
            else
                rep = zeros(c+1,c+1,r-1);
                for i=1:r-1
                    [rot,trans] = Tools.Trans_Mat(ob1.points(i,:),ob1.points(i+1,:));
                    rep(:,:,i)  = [rot,trans;...
                        zeros(1,c),1];
                end
            end
            ob1.reps = rep;
        end
        
        
        
        function ob1 = inv_rep(ob1)
            % Description: reconstructs points from representation and starting point
            % INPUT      : ob1 - instance of curve object
            % OUTPUT     : ob1 - a copy of the input object with reconstructed/updated "points" property
            
            
            if isempty(ob1.reps)
                error('There is no representation to invert !!!');
            else
                if ob1.type == 1
                    [~,c,r] = size(ob1.reps);
                    l = r;
                else
                    [~,c,r] = size(ob1.reps);
                    l = r+1;
                end
                im = zeros(l,c);
                im(1,:) = [ob1.start_point,1];
                for j=2:l
                    im(j,:)=(ob1.reps(:,:,j-1)*im(j-1,:)')';
                end
                
                ob1.points = im(:,1:c-1);
            end
        end
        
        
        
        function  ob2 = correspondence(ob1,ob2)
            % Description: Curve correspondence estimation under uniform
            %              sampling
            % INPUT      : ob1 and ob2  -  two instance of Curve classes
            % OUTPUT     : ob2          -  a copy of the second argument with rearranged
            %                              points to match the first argument points.
            
            
            % precaution
            if ~isempty(ob1.reps) && ~isempty(ob2.reps)
                ob1 = inv_rep(ob1);
                ob2 = inv_rep(ob2);
            else
                if ~isempty(ob1.points) && ~isempty(ob2.points)
                    ob1 = curve_rep(ob1);
                    ob2 = curve_rep(ob2);
                else
                    error('Argument objects are not initalized');
                end
            end
            
            r = size(ob1.reps,3);
            % for closed curves
            if ob2.type == 1
                val = (1:r)';
                dist = zeros(r,2);
                temp1 = ob2.reps;
                temp2 = Lie.inver(temp1);
                
                for k=1:2
                    % test matching points for different ordering directions
                    if k == 2
                        temp1 = temp2;
                        val = (1:r)';
                    end
                    for j=1:r
                        % 1-cyclic permutation
                        val=circshift(val,1);
                        trial = Lie.arrange(temp1,val);
                        dist(j,k) = Lie.dist_pro_lie(ob1.reps,trial);
                    end
                end
                [val2,x]= min(dist);
                [~,y]=min(val2);
                
                % rearranging the points in inverted ordering direction
                if y==2
                    ind = [x(1,2),2];
                    val = circshift((1:r)',ind(1,1));
                    n = flipud(ob2.points);
                    ob2.points = n(val,:);
                    % rearranging the points in the same ordering direction
                else
                    ind = [x(1,1),1];
                    val = circshift((1:r)',ind(1,1));
                    ob2.points = ob2.points(val,:);
                end
                % for open curves
            else
                dist = zeros(2,1);
                temp = Lie.inver(ob2.reps);
                
                dist(1,1) = Lie.dist_pro_lie(ob1.reps,ob2.reps);
                dist(2,1) = Lie.dist_pro_lie(ob1.reps,temp);
                
                if dist(1,1) > dist(2,1)
                    ob2.points = flipud(ob2.points);
                end
            end
            
            ob2 = curve_rep(ob2);
            
        end
        
        
        
        function ob2 = align_curves(ob1,ob2)
            % Description: Rotationally aligns the second argument with the
            %              first.
            %              WARNING: The function assumes point correspondance.
            %              -------
            %
            % INPUT      : ob1 and ob2 - two instances of Curve class
            % OUTPUT     : ob2         - rotationally aligned copy of the second argument.
            
            try
                if isempty(ob1.points) || isempty(ob2.points)
                    ob1 = curve_rep(ob1);
                    ob2 = curve_rep(ob2);
                end
            catch
                error('points property is not initalized');
            end
            
            % rotating around the fixed starting point
            centered_data1 = ob1.points - repmat(ob1.start_point,length(ob1.points),1);
            centered_data2 = ob2.points - repmat(ob2.start_point,length(ob2.points),1);
            
            rot = Tools.Opt_rot(centered_data2, centered_data1);
            
            centered_data2 = (rot*centered_data2')';
            ob2.points     =  centered_data2 + repmat(ob2.start_point,length(ob2.points),1);
            
            ob2 = curve_rep(ob2);
        end
        
        
        
        function ob2 = general_align(ob1,ob2)
            % Description: Computes the best general alignment of ob2 to
            %              ob1 through optimal distance estimation between orbits of the Abelian group.
            %
            % INPUT      : ob1 and ob2  - two instances of Curve class
            % OUTPUT     : ob2          - a copy of ob2 aligned upto reflection and rotation
            %                             based on uniform sampling correspondance estimation.
            %
            
            
            l = length(ob1.AbelianR);
            orbit_tests = cell(l,1);
            dist = zeros(l,1);
            
            for i=1:l
                orbit_tests{i,1} = ob2;
                orbit_tests{i,1}.points = Tools.normal_points((ob1.AbelianR(:,:,i)*orbit_tests{i,1}.points')');
                orbit_tests{i,1} = curve_rep(orbit_tests{i,1});
                orbit_tests{i,1} = correspondence(ob1,orbit_tests{i,1});
                orbit_tests{i,1} = align_curves(ob1,orbit_tests{i,1});
                dist(i,1) = ob1 - orbit_tests{i,1};
            end
            
            [~,I] = sortrows(dist);
            ob2 = orbit_tests{I(1,1),1};
        end
        
        
        function [dist,pair_dist] = minus(ob1,ob2)
            % Description: Computes the distance between two curve
            %              representations.
            % INPUT      : ob1 and ob2  - two instances of Curve class
            % OUTPUT     : dist         - scalar distance between the two curves
            %            : pair_dist    - element-wise distance between the two curves
            
            if ~isequal(ob1.start_point,ob2.start_point)
                warning('The starting points of the two curves are not the same!!');
            end
            
            r = size(ob1.reps,3);
            dist =  0; pair_dist = ones(r,1);
            for j=1:r
                pair_dist(j,1) = (Lie.geo_dist_SE(ob1.reps(:,:,j),ob2.reps(:,:,j)));
                dist = dist + pair_dist(j,1)^2;
            end
            dist = (dist)^0.5;
        end
        
        
        
        function curves = geodesic_path(ob1, ob2, steps)
            % Description: Computes instances along the geodesic curve
            %              connecting ob1 curve to ob2 curve.
            % INPUT      : ob1 and ob2   - two instances of Curve class.
            %            : steps         - (s) scalar, it is number of stops along
            %                              the geodeisc curve.
            % OUTPUT     : curves        - ((s+1)x1) cell of Curve instances.
            
            
            if isempty(ob1.reps) || isempty(ob2.reps)
                ob1 = curve_rep(ob1); ob2 = curve_rep(ob2);
            end
            
            rep1 = ob1.reps;  rep2 = ob2.reps;
            [~,c,r] = size(rep2);
            
            % step rate along the geodesic curve
            timestep = 1/steps;
            
            step = 0;
            
            curves{1,1} = ob1;
            
            for k=1:steps
                rep_temp = [];
                step = step + timestep;
                for i=1:r
                    rot = real(rep1(1:c-1,1:c-1,i)*(rep1(1:c-1,1:c-1,i)'*rep2(1:c-1,1:c-1,i))^(step));
                    tra = rep1(1:c-1,c,i) + (rep2(1:c-1,c,i) - rep1(1:c-1,c,i))*(step);
                    
                    mat = [rot,tra;...
                        zeros(1,c-1),1];
                    rep_temp = cat(3,rep_temp,mat);
                end
                curves{k+1,1} = Curve([],ob1.type,ob1.start_point);
                curves{k+1,1}.reps = rep_temp;
                curves{k+1,1} = inv_rep(curves{k+1,1});
            end
            curves{steps+1,1} = ob2;
        end
        
        
        
        function result = apply_deformation(ob1,def)
            % Description: computes the action of 'def' on the curve
            %              object 'ob1' from the left.
            % INPUT      : ob1    - instances of Curve class.
            %            : def    - 3-D array of rigid-transformation matrices.
            % OUTPUT     : result - instance of Curve class.
            try
                if isempty(ob1.reps)
                    ob1 = curve_rep(ob1);
                end
            catch
                error('Initalize the points property of the curve object argument');
            end
            
            if iscell(def)
                t = size(def,1);
                result = cell(t,1);
                for i=1:t
                    cur = Curve([],ob1.type,ob1.start_point);
                    % deformation acting from the left
                    cur.reps = Lie.group_pro(def{t,1},ob1.reps);
                    
                    % UNCOMMET the following to test deformation acting from the right
                    %cur.reps = Lie.group_pro(ob1.reps,def{t,1});
                    
                    result{i,1} = inv_rep(cur);
                end
            else
                cur = Curve([],ob1.type,ob1.start_point);
                % deformation acting from the left
                cur.reps = Lie.group_pro(def,ob1.reps);
                
                % UNCOMMENT the following to test deformation acting from the right
                %cur.reps = Lie.group_pro(ob1.reps,def);
                result = inv_rep(cur);
            end
            
            
        end
        
        
        
        function def = extract_deformation(ob1,ob2)
            % Description: computes the deformation 'def' that acts on the curve
            %              object 'ob1' from the left to give curve 'ob2'.
            % INPUT      : ob1, and ob2  - two instances of Curve class
            % OUTPUT     : def           - 3-D array of rigid-transformation matrices.
            
            if isempty(ob1.reps) || isempty(ob2.reps)
                curve_rep(ob1);
                curve_rep(ob2);
            end
            
            rep1 = ob1.reps;
            rep2 = ob2.reps;
            
            r = size(rep1,3);
            def = [];
            for j=1:r
                % deformation that acts from the left
                temp = rep2(:,:,j)/rep1(:,:,j);
                
                % UNCOMMENT the following to test deformation acting from the right
                %temp = rep1(:,:,j)\rep2(:,:,j);
                
                def = cat(3,def,temp);
            end
            
        end
        
    end
    
    
    methods(Static)
        
        
        function mean = shape_mean(curves,varargin)
            % Description: computes mean representation for a given set of
            %              Curve objects a dataset.
            %
            % INPUT      : curves      - (nx1) cell of Curve objects
            %            : [optional]  - (nx1) weight vector to compute a weighted mean.
            % OUTPUT     : mean        - karcher mean of the "curves"; instance of Curve object.
            
            
            [r,~] = size(curves);
            if nargin > 1
                weight = varargin{1};
            else
                weight = repmat(1/r,r,1);
            end
            
            
            R = curves{1,1}.reps;
            [~,~,n] = size(R);
            
            M_rep = [];
            for j=1:n
                Temp = [];
                for i=1:r
                    Temp= cat(3,Temp,curves{i,1}.reps(:,:,j));
                end
                M_rep = cat(3, M_rep, Lie.meanSE(Temp,weight));
            end
            
            mean =  Curve([], curves{1,1}.type,curves{1,1}.start_point);
            mean.reps = M_rep;
            mean = inv_rep(mean);
        end
        
        
        
        function [means,ind] = kmeans_shape(curves,K,itt,means)
            % Description: computes K-clusters of a curved shape dataset.
            % INPUT      : curves - (nx1) cell of Curve objects
            %            : K      - number of clusters
            %            : itt    - number of maximum itterations.
            %            : means  - (Kx1) cell of Curve objects; inital
            %                       means.
            % OUTPUT     : means  - (Kx1) cell of Curve objects; estimated means.
            %            : ind    - (Kx1) cell of vectors specifying the index
            %                        of the each clsuters' elements.
            
            
            [r,~] = size(curves);
            
            if r == 0
                means = []; ind = []; return;
            end
            
            % select inital mean shape if not given
            for k=1:K
                if isempty(means{k,1})
                    h = randi(r,1);
                    means{k} = curves{h,1};
                end
            end
            
            dist = zeros(K,r);
            
            
            tempIn = []; tempmeans = [];
            
            objective_cost = zeros(itt,1);
            %-- EM algorithm
            for q=1:itt
                flag = 0;
                
                %E-step
                for i=1:r
                    for k=1:K
                        dist(k,i) = means{k} - curves{i,1};
                    end
                end
                [~,In] = min(dist);
                
                %M-step
                for j=1:K
                    d1 = sum(In(:)==j);
                    val = find(In==j);
                    % each mean element has to have at least 3 elements
                    if d1 <= 2
                        flag=1;
                        break;
                    end
                    means{j} = Curve.shape_mean(curves(val'),repmat(1/d1,d1,1));
                    
                    % objective cost
                    objective_cost(q,1) = objective_cost(q,1) + sum(dist(j,val).^2);
                end
                
                if flag
                    if isempty(tempIn)
                        % re-initalize and try again
                        for i=1:K
                            means{i} = curves{randi(r,1),1};
                        end
                        continue;
                    else
                        means = tempmeans;
                        In = tempIn;
                        break;
                    end
                end
                
                tempIn = In;
                tempmeans = means;
            end
            %-- end of EM
            
            % elemnts of each cluster
            for t=1:K
                ind{t,1} = find(In(:) == t);
            end
             
        end
        
        
        function plot_corr(ob1,ob2, varargin)
            % Description: plots point correspondance between curves
            %              ob1 and ob2
            % INPUT      : ob1, and ob2 - Two instances of Curve class
            
            
            c = size(ob1.points,2);
            if ob1.type == 0
                sh1 =  ob1.points;
                sh2 =  ob2.points;
            else
                sh1 =  [ob1.points;ob1.points(1,:)];
                sh2 =  [ob2.points;ob2.points(1,:)];
            end
            temp = [sh1(:)';sh2(:)'];
            temp = temp(:);
            
            
            Cor2 = reshape(temp,[2*length(sh1),c]);
            
            figure;
            if c == 3
                scatter3(sh1(:,1),sh1(:,2),sh1(:,3),'*b'); hold on;
                line(sh1(:,1),sh1(:,2),sh1(:,3),'color','blue'); hold on;
                scatter3(sh2(:,1),sh2(:,2),sh2(:,3),'*r'); hold on;
                line(sh2(:,1),sh2(:,2),sh2(:,3),'color','red');hold on;
                for i=1:2:(length(Cor2)-1)
                    line(Cor2(i:i+1,1),Cor2(i:i+1,2),Cor2(i:i+1,3),'LineStyle','--'); hold on;
                end
            else
                scatter(sh1(:,1),sh1(:,2),'*b'); hold on;
                line(sh1(:,1),sh1(:,2),'color','blue'); hold on;
                scatter(sh2(:,1),sh2(:,2),'*r'); hold on;
                line(sh2(:,1),sh2(:,2),'color','red'); hold on;
                for i=1:2:(length(Cor2)-1)
                    line(Cor2(i:i+1,1),Cor2(i:i+1,2),'LineStyle','--','color','black'); hold on;
                end
                
            end
            if nargin > 2
                if ischar(varargin{1})
                    title(varargin{1});
                end
            else
                title('Point correspondence estimation');
            end
            axis on; axis square; grid on; hold off;
        end
        
        
        
        function plot_curves(ob1,varargin)
            % Description: visualize a given Curve object/cell of curve
            %              objects.
            % INPUT      : ob1        - curve object or cell of curve objects.
            %            : [optional] - text as a title.
            
            [t,~] = size(ob1);
            s = 1;
            
            
            f1 = figure;
            
            for j=1:t
                if iscell(ob1)
                    set(f1,'OuterPosition',[200,200,1200,300]);
                    subplot(s,t,j);
                    Curve.plotme(ob1{j});
                else
                    Curve.plotme(ob1);
                    break;
                end
                
            end
            if nargin > 1
                if ischar(varargin{1})
                    suptitle(varargin{1});
                    return;
                else
                    warning('Title is not text, default title will be displayed.');
                end
            end
            
            if t > 1
                suptitle('Geodesic deformation from the first curve to the last');
            end
            
            hold off;
        end
        
        
        
        function plotme(ob1)
            % Description: visualize a given Curve object
            % INPUT      : object of curve class
            
            try
                if isempty(ob1.points)
                    ob1 = inv_rep(ob1);
                end
            catch
                error('First initalize points OR reps property of the argument.');
            end
            
            points_t = ob1.points;
            
            % scale point for clear visualization
            %dec = 0;
            %points_t = round(points_t.*10^(dec))./10^(dec-1);
            
            r = size(points_t,2);
            
            if ob1.type == 1
                points_t = [points_t;points_t(1,:)];
            end
            
            if r == 2
                % starting point
                scatter(points_t(1,1),points_t(1,2),'filled','g'); hold on;
                % end point
                scatter(points_t(end-3,1),points_t(end-3,2),'filled','r'); hold on;
                line(points_t(:,1),points_t(:,2),'LineWidth',2,'color','black');
                axis off;
            else
                % starting point
                scatter3(points_t(1,1),points_t(1,2),points_t(1,3),'filled','g'); hold on;
                % end point
                scatter3(points_t(end,1),points_t(end,2),points_t(end-1,3),'filled','r'); hold on;
                line(points_t(1:end,1),points_t(1:end,2),points_t(1:end,3),'LineWidth',2, 'color','black');
                % view([35 -19]);
                view([-37,32]);
                axis equal; grid off;
                box on;
            end
            
        end
        
    end
end
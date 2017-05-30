

% Description: DP_OPTIMAL_POINT_SAMPLING_SINGLE samples c2 optimally.
%
%
% INPUT      : c1     - uniformly sampled instance of a Curve class.
%            : c2     - instance of a curve with more points than c1, i.e.
%                       |c2.points| >> |c1.points|
%            : alpha  - a scalar constant to weight deformation term.
%            : beta   - a scalar constant to weight constraint term.
%            : win    - a positive scalar determing the window size of the
%                        charts.
% OUTPUT     : index3 - index of the choosen points from c2.points.

% HOWTO      : see Demo_HOWTO_curve_representation.m and Tools.DP_sampling

% 2017  Girum G. demisse, girumdemisse@gmail.com/girum.demisse@uni.lu
%       Computer vision team, University of Luxembourg.
%------------------------------------------------


function index3 = dp_optimal_point_sampling_Single(c1,c2,alpha,beta,win)
    
    
    % defines the radius of the window size
    if win <= 0
        error('Window size cannot be negative !!');
    end
    
    [r,~] = size(c1.points);
    [c,~] = size(c2.points);
    
    
    % Variables for intermidate computations
    cost1 = Inf(c,r);
    cost2 = Inf(c,r);
    cost3 = Inf(c,r);
    
    temp_ind1 = zeros(c,r);
    temp_ind2 = zeros(c,r);
    temp_ind3 = zeros(c,r);
    
    cost1(1,1) = 0;
    cost2(1,1) = 0;
    cost3(1,1) = 0;
    
    % Three solution paths
    index3 =  zeros(r,1);
    index2 = zeros(r,1);
    index1 = zeros(r,1);
    
    
    % Initalization
    index1(1,1) = 1;
    index2(1,1) = 1;
    index3(1,1) = 1;
    
    temp_ind1(1,1) = 1;
    temp_ind2(1,1) = 1;
    temp_ind3(1,1) = 1;
    
    
    sliding_rate = ceil(c/(r));
    
    % initalize constraint term
    if size(c1.points,2) == 2
         % area based constraint for planar curved shapes 
         [~,uni_rate]  = Tools.Curve_area(c2.points(1:sliding_rate:c,:),c2.type);
    else
         % length based constraint for curved shapes in > 2 dimensional space
         [~,uni_rate]  = Tools.Curve_length(c2.points(1:sliding_rate:c,:),c2.type);
    end
  
    % solve for each sub-functions
    for i=1:(r-1)
        p = sliding_rate*(i-1);
        start = max(2,p-win);
        las = min(c,win+p);
        
        pre = max(2,start-sliding_rate);
       
        if i == 1
            pre = 1;
        end
        % solving for the minimum and minimizer of the (cost matrix)
        for j=start:las  
              
              [temp1,temp2,temp3] = opt_subStructure(j,c2,c1,pre,uni_rate,alpha,beta,i);
              [cost1(j,i+1),temp_ind1(j,i+1)] = min(cost1(:,i) + temp1);
              [cost2(j,i+1),temp_ind2(j,i+1)] = min(cost2(:,i) + temp2);
              [cost3(j,i+1),temp_ind3(j,i+1)] = min(cost3(:,i) + temp3);  
           
        end
    end
    
    
    % Work backwards to get the optimal sampling solution
    [~,index1(r,1)] = min(cost1(:,r));
    [~,index2(r,1)] = min(cost2(:,r));
    [~,index3(r,1)] = min(cost3(:,r));
    
    for i=r:-1:2
        % cost of deformation
        index1(i-1,1) = temp_ind1(index1(i,1),i);
        % cost of the area
        index2(i-1,1) = temp_ind2(index2(i,1),i);
        % cost of constrained objective functional
        index3(i-1,1) = temp_ind3(index3(i,1),i);
        
    end
    
    % Ploting the cost matrix and the three solutions
     Tools.plot_opt_paths(cost3,index1,index2,index3);
end


function [cost1,cost2,cost3] = opt_subStructure(j,c2,c1,s,uni_rate,alpha,beta,F)
    
    len = size(c2.points,1);
    cost1 = inf(len,1);
    cost2 = inf(len,1);
    cost3 = inf(len,1);
    fix = c1.reps(:,:,F);
    
    for i=s:(j-1)
        [rot,trans] = Tools.Trans_Mat(c2.points(i,:),c2.points(j,:));
        MAT =  [rot,trans;...
            zeros(1,length(c2.start_point)),1];
        
        %**** FIXME: avoid checking this condition in each itteration
        if size(c2.points,2) == 2 
            % area based constraint for planar curved shapes 
            
            dx = c2.points(i,1)-c2.points(j,1);
            dy = 0.5*(c2.points(i,2)+c2.points(j,2));

            % NOTE the absolute value is used to ensure that this is a cost and
            % and the lowest possible value is zero.

            constraint = abs(abs(uni_rate(F)) - abs(dx*dy));
            
                % closing constraint
                if (F == (length(uni_rate)-1)) && (c2.type == 1) 
                    dx = c2.points(j,1)-c2.points(1,1);
                    dy = 0.5*(c2.points(j,2)+c2.points(1,2));

                    constraint = constraint + abs(abs(uni_rate(F+1)) - abs(dx*dy));
                end
        else
            % length based constraint for curved shapes in > 2 dimensional space
            
            l_segment = norm(c2.points(i,1)-c2.points(j,1));
            constraint = abs(abs(uni_rate(F)) - l_segment);
        end
        
        cost1(i,1) =  alpha*(Lie.geo_dist_SE(fix,MAT))^2;
        cost2(i,1) =  beta*constraint;
        
                 % closing deformation cost
                if (F == (length(uni_rate)-1)) && (c2.type == 1)
                     [rot,trans] = Tools.Trans_Mat(c2.points(j,:),c2.points(1,:));
                     MAT2 =  [rot,trans;...
                     zeros(1,length(c2.start_point)),1]; 
                   
                     cost1(i,1) = cost1(i,1) + alpha*(Lie.geo_dist_SE(c1.reps(:,:,F+1),MAT2))^2;
                end
        
        
        % constrained cost functional
        cost3(i,1) =  cost1(i,1) + cost2(i,1);
    end 
end

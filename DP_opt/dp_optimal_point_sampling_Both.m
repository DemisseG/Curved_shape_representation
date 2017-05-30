
% Description: DP_OPTIMAL_POINT_SAMPLING_SINGLE samples c1 and c2 optimally.
%
%
% INPUT      : c1     - instance of a Curve class.
%            : c2     - instance of a Curve class, NOTE: |c2.points| == |c1.points|
%            : alpha  - a scalar constant to weight deformation term.
%            : beta   - a scalar constant to weight constraint term.
%            : win    - a positive scalar determing the window size of the
%                        charts.
%            : sample - scalar value, determins the number of points to be sampled
%                       sample << |c2.points|.

% OUTPUT     : index3 - (sample,2) index of the choosen points from
%              c1.points, index3(:,1), and from c2.points, index3(:,2).

% HOWTO      : see Demo_HOWTO_curve_representation.m and Tools.DP_sampling

% 2017  Girum G. demisse, girumdemisse@gmail.com/girum.demisse@uni.lu
%       Computer vision team, University of Luxembourg.
%------------------------------------------------


function index3 = dp_optimal_point_sampling_Both(c1,c2,alpha,beta,win,sample)
    
    
    if win <= 0
        error('Window size cannot be negative !!');
    end
    
    
    [r,~] = size(c1.points);
    [c,~] = size(c2.points);
    
    
    
    % Variables for intermidate computations
    cost1 = Inf(r,c,sample);
    cost2 = Inf(r,c,sample);
    cost3 = Inf(r,c,sample);
    
    temp_ind1 = zeros(r,c,2,sample);
    temp_ind2 = zeros(r,c,2,sample);
    temp_ind3 = zeros(r,c,2,sample);
    
    cost1(1,1,1) = 0;
    cost2(1,1,1) = 0;
    cost3(1,1,1) = 0;
    
    % Three solution paths
    index1 = zeros(sample,2);
    index2 = zeros(sample,2);
    index3 =  zeros(sample,2);
    
    % Initalization (the first points match)
    index1(1,1:2) = [1,1];
    index2(1,1:2) = [1,1];
    index3(1,1:2) = [1,1];
    
    temp_ind1(1,1,:,1) = [1,1];
    temp_ind2(1,1,:,1) = [1,1];
    temp_ind3(1,1,:,1) = [1,1];
    
    
    sliding_rate = ceil(c/sample);
    
    % area constraint: SHOULD NOT BE Normalized.
    if size(c1.points,2) == 2
        % area based constraint for planar curved shapes
        [~,uni_rate1]  = Tools.Curve_area(c1.points(1:sliding_rate:c,:),c1.type);
        [~,uni_rate2]  = Tools.Curve_area(c2.points(1:sliding_rate:c,:),c2.type);
    else
        % length based constraint for curved shapes in > 2 dimensional space
        [~,uni_rate1]  = Tools.Curve_length(c1.points(1:sliding_rate:c,:),c1.type);
        [~,uni_rate2]  = Tools.Curve_length(c2.points(1:sliding_rate:c,:),c2.type);
    end
    
    % solve for each sub-functions
    for i=1:(sample-1)
        
        % sliding window of search size
        p = sliding_rate*(i-1);
        start = max(2,p-win);
        las = min(c,win+p);
        
        pre = max(2,start-sliding_rate);
        
        if i==1
            pre = 1;
        end
        % iterating through the first dimension (first curve points)
        for j1 =start:las          
            % iterating through the second dimension(second curve points)
            for j2 = start:las
                [temp1,temp2,temp3] = opt_subStructure(j1,j2,c2,c1,pre,uni_rate1,uni_rate2,alpha,beta,i);

                [cost1(j1,j2,i+1),temp_ind1(j1,j2,:,i+1)] = Tools.Min_M(cost1(:,:,i) +temp1);
                [cost2(j1,j2,i+1),temp_ind2(j1,j2,:,i+1)] = Tools.Min_M(cost2(:,:,i) +temp2);
                [cost3(j1,j2,i+1),temp_ind3(j1,j2,:,i+1)] = Tools.Min_M(cost3(:,:,i) +temp3);
            end
        end
    end
    
    
    % Work backwards to get the optimal sampling solution
    [~,index1(sample,:)] = Tools.Min_M(cost1(:,:,sample));
    [~,index2(sample,:)] = Tools.Min_M(cost2(:,:,sample));
    [~,index3(sample,:)] = Tools.Min_M(cost3(:,:,sample));
    for i=sample:-1:2
        % cost of deformation
        rr =index1(i,1); cc=index1(i,2);
        index1(i-1,:) = temp_ind1(rr,cc,:,i);
        
        % cost of the area
        rr2=index2(i,1); cc2=index2(i,2);
        index2(i-1,:) = temp_ind2(rr2,cc2,:,i);
        
        % cost of the conjugate effect
        rr3=index3(i,1); cc3=index3(i,2);
        index3(i-1,:) = temp_ind3(rr3,cc3,:,i);
    end
    
    % Ploting the cost matrix and the three solutions
    Tools.plot_opt_paths(cost3,index1,index2,index3);
   
end



function [cost1,cost2,cost3] = opt_subStructure(j,j2,c2,c1,s,uni_rate1,uni_rate2,alpha,beta,F)
    
    r = length(c1.points); c = length(c2.points);
    cost1 = inf(r,c);
    cost2 = inf(r,c);
    cost3 = inf(r,c);
    
    
    d = size(c2.start_point,2);
    for i=s:(j-1)
        
        %**** FIXME: avoid checking this condition in each itteration
        
        if size(c1.points,2) == 2
            % area based constraint for planar curved shapes
            dx = c1.points(i,1)-c1.points(j,1);
            dy = 0.5*(c1.points(i,2)+c1.points(j,2));

            constraint_1 = abs(abs(uni_rate1(F)) - abs(dx*dy));
                % closing constraint of area
                if (F == (length(uni_rate1)-1)) && (c1.type == 1) 
                    dx = c1.points(j,1)- c1.points(1,1);
                    dy = 0.5*(c1.points(j,2)+c1.points(1,2));

                    constraint_1 = constraint_1 + abs(abs(uni_rate1(F+1)) - abs(dx*dy));
                end
                % end of closing constarint
        else
            % length based constraint for curved shapes in > 2 dimensional space
            l_segment = norm(c1.points(i,1)-c1.points(j,1));
            constraint_1 = abs(abs(uni_rate1(F)) - l_segment);
        end
        
        [rot,trans] = Tools.Trans_Mat(c1.points(i,:),c1.points(j,:));
        MAT1 =  [rot,trans;...
                zeros(1,d),1];
           
         % closing deformation constraint for first curve
         if (F == (length(uni_rate1)-1)) && (c1.type == 1)
             [rot,trans] = Tools.Trans_Mat(c1.points(j,:),c1.points(1,:));
             MAT1_C =  [rot,trans;...
             zeros(1,length(c1.start_point)),1]; 
         else
             MAT1_C = eye(length(c1.start_point)+1);
         end
         
        for k=s:(j2-1)
            
            if size(c2.points,2) == 2
                dx2 = c2.points(k,1)-c2.points(j2,1);
                dy2 = 0.5*(c2.points(k,2)+c2.points(j2,2));

                constraint_2 = abs(abs(uni_rate2(F)) - abs(dx2*dy2));
                % closing constraint of area
                if (F == (length(uni_rate2)-1)) && (c2.type == 1) 
                    dx = c2.points(j2,1)-c2.points(1,1);
                    dy = 0.5*(c2.points(j2,2)+c2.points(1,2));

                    constraint_2 = constraint_2 + abs(abs(uni_rate2(F+1)) - abs(dx*dy));
                end
                % end of closing constraint
            else
                l_segment = norm(c2.points(k,1)-c2.points(j2,1));
                constraint_2 = abs(abs(uni_rate2(F)) - l_segment);
            end
            
            [rot2,trans2] = Tools.Trans_Mat(c2.points(k,:),c2.points(j2,:));
            
            MAT2 =  [rot2,trans2;...
                     zeros(1,d),1];
            
            % closing deformation constraint for second curve
             if (F == (length(uni_rate2)-1)) && (c2.type == 1)
                 [rot,trans] = Tools.Trans_Mat(c2.points(j2,:),c2.points(1,:));
                 MAT2_C =  [rot,trans;...
                 zeros(1,length(c2.start_point)),1]; 
             else
                 MAT2_C = eye(length(c2.start_point)+1);
             end
                 
                 
                 
                 
            % cost of deformation
            cost1(i,k) = alpha*((Lie.dist_pro_lie(MAT1,MAT2))^2) + alpha*((Lie.dist_pro_lie(MAT1_C,MAT2_C))^2);
            % cost of area constraint
            cost2(i,k) =  beta*(constraint_1 + constraint_2);
            % constrained cost functional
            cost3(i,k) =  cost1(i,k) + cost2(i,k);
        end
    end
end

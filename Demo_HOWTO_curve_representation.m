

% DEMO script on how to use the Curve class  
%          * to build curve representation.
%          * to compute point correspondance (using both uniform and
%            optimal sampling).
%          * to compute geodesic distance between curves.
%          * to plot results.
%-------------------------------------------------
% HINT:   Simple approach to tune the free parameters effeciently
%          1) assign large window_size(>10 or better >= 20).
%          2) fix alpha to 1 and tune beta in [0,1].
%          3) In case of optimally sampling both shapes,
%             it is better to tune alpah and beta together.
%
%         *** To achieve high accuracy in retrival increse resolution, that is, 
%             set: num_points(>=200) and sample_points(>=100), and set
%             window_size >= 20. However, this entails high computational time. ***
%--------------------------------------------------     
% NOTICE:  Implementation specific notes:
%          1) mod(num_points,sample_points) should be 0.
%          2) The starting points are aligned using
%             uniform sampler for optimal sampling, for high accuracy use
%             optimal sampling to estimate starting point.

%--------------------------------------------------
% 2017  Girum G. Demisse, girumdemisse@gmail.com/girum.demisse@uni.lu
%       Computer vision team, University of Luxembourg.
%--------------------------------------------------



clear;

% parameters
num_points = 200;              % number of points to extract from data source. 
sample_points = 100;           % number of points to be sampled (either uniformly or optimally).  
window_size = 10;              % search space size to sample a point
alpha = 1;                     % weighting parameter for the deformation cost
beta = 0.1;                    % weighting parameter for the area cost
D_type  = 1;                   % = 1 fix the first argument and optimally sample the second argument.
                               % = 0 optimally sample both arguments, for quick
                               %     test, set window_size=5.
  
  

%-- load two ordered data points
D = load('KIMIA99'); D = D.KIMIA99;
%D = load('KIMIA216'); D = D.KIMIA216;

%-- select points
points1 = D{2,7}; points2 = D{3,7};

%-- sample and normalize points
points1 = Tools.process(points1,num_points);
points2 = Tools.process(points2,num_points);


%-- initalize Curve class object
c1 = Curve(points1); 
c2 = Curve(points2);

%-------------------------------------------------------------------------------------------------------------
%-- General alignment ESTIMATION of c2 to c1 is based on point correspondance result from uniform sampling.
%-- CAUTION: the alignment result might not be optimal (specially in large deformations) since the 
%            point correspondance is based on uniform sampling.            
%-------------------------------------------------------------------------------------------------------------
if size(c2.points,2) == 2
   % for planar curves
   c2 = general_align(c1,c2);
else
   % for non-planar curves
   c2 = correspondence(c1,c2); 
   c2 = align_curves(c1,c2);
end


%--------------------------------------------------------------
% --- UNCOMMENT the following to test the symmetric nature of sampling both curves
% temp = c2;
% c2 = c1;
% c1 = temp;
%---------------------------------------------------------------

%-- Correspondance based on optimal sampling of points
[c1_r,c2_r] = Tools.DP_sampling(c1,c2,alpha,beta,window_size,D_type,sample_points);           

%--- Correspondance based on uniform sampling of points
c1_U = c1; 
c1_U.points = Tools.process(c1.points,sample_points); 
c1_U = curve_rep(c1_U);

c2_U = c2; 
c2_U.points = Tools.process(c2.points,sample_points); 
c2_U = curve_rep(c2_U);

%--- Geodesic curve between the curves with optimal sampling 
curves = geodesic_path(c1_r,c2_r,5);
Curve.plot_curves(curves,'Based on optimal sampling');

%--- Geodesic curve between the curves with uniform sampling 
curves2 = geodesic_path(c1_U,c2_U,5);
Curve.plot_curves(curves2,'Based on uniform sampling');

%--- Plot estimated point correspondance 
Curve.plot_corr(c1_r,c2_r,'Optimal sampling based');
Curve.plot_corr(c1_U,c2_U,'Uniform sampling based');

%--- Geodesic distance
dist_opt = c1_r - c2_r;
disp('Geodesic distance under optimal sampling:'); disp(dist_opt);
dist_uni = c1_U - c2_U;
disp('Geodesic distance under uniform sampling:'); disp(dist_uni);

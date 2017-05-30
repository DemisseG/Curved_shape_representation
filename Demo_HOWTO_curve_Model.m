
 
% DEMO script on How to estimate means and/or clusters of shapes. 
% For further detail Consult: "Deformation Based Curved shape Representation", submitted to TPAMI.
%                              and Demo_HOWTO_curve_representation.m for detail description of the
%                              parameters.
%
% 2017  Girum G. Demisse, girumdemisse@gmail.com/girum.demisse@uni.lu
%       Computer vision team, University of Luxembourg.
%------------------------------------------------


%clear;

%-- parameters
num_points = 100;
sample_points = 50;
window_size = 10; 
alpha = 1; 
beta = 0.3; % use to 0.3

%-- load dataset
data =  load('KIMIA99');
D = data.KIMIA99;

%-- curved shape category to be modelled
class_type = 1;

%-- training data size
last = size(D,1);

%-- select uniformly sampled template curved shape to align all the other dataset.
points = Tools.process(D{4,class_type},num_points);
temp = Curve(points);

%-- uniform sampling and normalization of temp to temp_U

temp_U = temp;
l = ceil(num_points/sample_points);
temp_U.points = Tools.normal_points(temp_U.points(1:l:end,:));
temp_U = curve_rep(temp_U);
 
dynamic = 1; %  set to 1  to enable dynamic programming for optimal sampling
curves = cell(last,1);

for i=1:last
    
    points_t = Tools.process(D{i,class_type},num_points);
    cur = Curve(points_t);
    
    %-- General alignment ESTIMATION of c2 to c1. 
    cur = general_align(temp,cur); 
 
    %--- UNCOMMENT to use rotation only based alignment
    % cur = correspondence(temp,cur);
    % cur = align_curves(temp,cur);
    
    if dynamic   
       opt_index  = dp_optimal_point_sampling_Single(temp_U,cur,alpha,beta,window_size);
       cur.points = Tools.normal_points(cur.points(opt_index,:));
       cur =  curve_rep(cur);
    end
    % collect aligned and corresponding curved shapes  
    curves{i,1} =  cur;
end


%-- UNCOMMENT the following to compute mean of curves dataset
% mean_s = Curve.shape_mean(curves,repmat(1/length(curves),length(curves),1));
% Curve.plot_curves(mean_s);

%-- clustering parameters
K = 2;    % number of clusters
itt = 10; % number of iterations

%-- inital guess of the mean observations
inital = {curves{2,1};curves{6,1}}; % second curve use to 12

%-- compute K-clusters based on EM
[means,elements] = Curve.kmeans_shape(curves,K,itt,inital);

%-- viaualize the estimated clusters with their means.
visualizeK_means; 




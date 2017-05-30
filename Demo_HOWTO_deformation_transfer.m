

% DEMO script on How to Transport deformation
% For Further detail Consult: "Deformation Based Curved shape Representation", submitted to TPAMI.
%                              and Demo_HOWTO_curve_representation.m details on the
%                              parameters.
%
% 2017  Girum G. Demisse, girumdemisse@gmail.com/girum.demisse@uni.lu
%       Computer vision team, University of Luxembourg.
%------------------------------------------------


clear;

%-- parameters
num_points = 100;
sample_points = 50;
window_size = 10; 
alpha = 1; 
beta = 0.1;
D_type = 1;

%-- load data
D = load('KIMIA99'); D = D.KIMIA99;

%-- first curve
points1 = D{1,7};
%-- first curve under deformation
points2 = D{10,7};
%-- second curve that is similar to the first but NOT identical
points3 = D{2,7};

%-- initalize curve class
points1 = Tools.process(points1,num_points);
points2 = Tools.process(points2,num_points);
points3 = Tools.process(points3,num_points);

c1 = Curve(points1);
c2 = Curve(points2);
c3 = Curve(points3);

%-- General alignment ESTIMATION of c2 to c1.
c2 = general_align(c1,c2);

%-- General alignment ESTIMATION of c3 to c1. 
c3 = general_align(c1,c3);

%-- Fix c1 to uniform sampling and sample c2 optimally.
[c1_U,c2_O] = Tools.DP_sampling(c1,c2,alpha,beta,window_size,D_type,sample_points);

%-- Fix c1 to uniform sampling and sample c3 optimally.
[~,c3_O] = Tools.DP_sampling(c1,c3,alpha,beta,window_size,D_type,sample_points);


%---------------------------------------
% Deformation Transport
%----------------------------------------

%-- extracting the deformation from c1_U to c3_O
def = extract_deformation(c1_U,c2_O);

%-- Left action of def on c2_O
c4 = apply_deformation(c3_O,def);


%-- Result plots

f1 = figure;
set(f1,'OuterPosition',[400,400,700,400]);
sh1 = subplot(2,2,1); title(sh1,'c1');
sh2 = subplot(2,2,2); title(sh2,'Observed deformation of c1');
sh3 = subplot(2,2,3); title(sh3,'c3, similar but not identical to c1');
sh4 = subplot(2,2,4); title(sh4,'Transfered deformation to c3');

Curve.plot_curves(c1_U,'c1');
fig = gcf;
copyobj(allchild(get(fig,'CurrentAxes')),sh1); 
close(fig);

Curve.plot_curves(c2_O,'Observed deformation of c1.')
fig = gcf;
copyobj(allchild(get(fig,'CurrentAxes')),sh2); 
close(fig);

Curve.plot_curves(c3_O,'c3, similar but not identical'); 
fig = gcf;
copyobj(allchild(get(fig,'CurrentAxes')),sh3); 
close(fig);

Curve.plot_curves(c4,'Transfered deformation to c3.'); 
fig = gcf;
copyobj(allchild(get(fig,'CurrentAxes')),sh4); 
close(fig);

%--- deformation alng the geodesic path and 
curves1 = geodesic_path(c1_U,c2_O,5);
Curve.plot_curves(curves1,'Geodesic deformation from c1 to c3.');

curves2 = geodesic_path(c3_O,c4,5);
Curve.plot_curves(curves2,'Geodesic deformation from c2 to its predicted deformation.');



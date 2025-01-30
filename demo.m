close all

%
% DEFINE constants/fxns here
%

n = 6;
c_v = 0.2; 
c_i = 0.2; 
s = 0.05; 

V = @(j) 1/(1+exp(-10*(j-3)/6));
ben_v = @(j, n) (V(j)-V(0))/(V(n)-V(0)); 
ben_i = @(j, n) (V(j)-V(0))/(V(n)-V(0)); 

x0 = [0.3,0.2,0.5];
tF = 50;

%
% Instantiate class
%
A = nonlinear_pure(n, c_v, c_i, s, ben_v, ben_i);

%
% Demo all features
%
black = [0,0,0]; red = [1,0,0]; blue = [0,1,0]; green = [0,0,1];

labels = {'A','B','C'};

%% plot total velocity magnitude, with quivers
figure(1); hold on;
isomatrix_velocity(A); % This takes a few seconds to run
colorbar;

add_gridlines(10); % 10 gridlines in each direction
isomatrix_fixedpointedge(A);
isomatrix_quiver(A);
add_labels(labels);

%% plot a trajectory on simplex starting at x0 until time=tF
isomatrix_trajectory(A,x0,tF,'Color',red);

%% line plot of same simulation
figure(2);
line_plot(A,x0,tF,'Labels',labels);
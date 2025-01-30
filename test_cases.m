% Linear Benefits, Archetype 1 clearly wins

n = 6; 
c_v = 0.1; 
c_i = 0.1;
s = 1; 
ben_v = @(j, n) j / n ; 
ben_i = @(j, n) j / n ; 

x0 = [0.3,0.2,0.5];
tF = 10;

% Nonlinear Benefits, stable A1 - A3 equilibrium and stable A2-A3
% equilibrium
n = 6;
c_v = 0.2; 
c_i = 0.2; 
s = 0.0; 

V = @(j) 1/(1+exp(-10*(j-3)/6));
ben_v = @(j, n) (V(j)-V(0))/(V(n)-V(0)); 
ben_i = @(j, n) (V(j)-V(0))/(V(n)-V(0)); 

x0 = [0.3,0.2,0.5];
tF = 50;

% All have equivalent fitness of 0
n = 4; 
c_v = 0;
c_i = 0; 
s = 0; 
ben_v = @(j, n) 0 ; 
ben_i = @(j, n) 0 ; 

x0 = [0.1,0.2,0.7];
tF = 10;


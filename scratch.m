close all

% Instantiate
n = 6;
c_v = 0.2; 
c_i = 0.2; 
s = 0.0; 

V = @(j) 1/(1+exp(-10*(j-3)/6));
ben_v = @(j, n) (V(j)-V(0))/(V(n)-V(0)); 
ben_i = @(j, n) (V(j)-V(0))/(V(n)-V(0)); 

x0 = [0.5,0,0.5];
tF = 10;

A = nonlinear_pure(n, c_v, c_i, s, ben_v, ben_i);

% Define the x range
x_values = linspace(0, 1, 100); % 100 points between 0 and 1

% Preallocate arrays for W1 and W2
y1 = zeros(size(x_values));
y2 = zeros(size(x_values));
y3 = zeros(size(x_values));

for i = 1:length(x_values)
    x = x_values(i); 
    triple = [x, 1-x, 0]; 
    y1(i) = A.W1(triple);
    y2(i) = A.W2(triple); 
    y3(i) = A.W3(triple); 
end

% Plot the functions
figure; 
plot(x_values, y1, 'Color', [0 0.4470 0.7410], 'LineWidth', 2); 
hold on; 
plot(x_values, y2, 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2); 
plot(x_values, y3, 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2); 
hold off; 

% Add labels, title, and legend
xlabel('x'); % Label for x-axis
ylabel('Fitness'); % Label for y-axis
title('A1, A2, A3 Fitness'); % Title of the graph
legend('A1', 'A2', 'A3', 'Location', 'best'); % Add legend
grid on; % Add a grid for better readability
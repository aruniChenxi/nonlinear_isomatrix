function [] = isomatrix_fixedpointedge(A, index, varargin)

p = inputParser;
effective_index = 1;
% if odd number of arguments:
if ((mod(nargin, 2) == 1) && (nargin > 2))
    % no user specified index, with extra arguments
    varargin = [{index}, varargin];
    index = 1;
    effective_index = 0;
    addOptional(p, 'index', index);
elseif (nargin == 1)
    % no user specified index, with no extra arguments
    index = 1;
    effective_index = 0;
end

%% set up default values for optional parameters: ('Color' and 'Labels')
color = [0, 0, 0];
labels = {'', '', ''};

% validation of user input color:
vectorValidator = @(x) validateattributes(x, {'numeric'}, {'size', [1, 3]});
addParameter(p, 'Color', color, vectorValidator)

% read in optional parameters
[nParams] = length(varargin);
for param = 1:1:(nParams / 2)
    ind = (param - 1) * 2 + 1;
    if strcmp(varargin{ind}, 'Color')
        color = varargin{ind+1};
    elseif strcmp(varargin{ind}, 'Labels')
        labels = varargin{ind+1};
        labelLength(labels);
        labelType(labels);
    end
end

h = gcf;
figure_number = h.Number;
figure(figure_number);
hold on;

%% fix margins if there are many games
H = (1 / 2) * tan(60*pi/180);
L = 1;
del = (1 - H) / 2;
mul = index;
xlim([-mul * del, L + mul * del]);
ylim([-mul * del, H + mul * del]);

% offset axes parameters (probably best not to change)
line = [0, 1];
REL_DIST = 0.05;
delta = [0, 0];

us_ms = 12;
s_ms = 50;

ylim([-REL_DIST * index, (1 + REL_DIST * index) / sin(pi/3)]);
xlim([-REL_DIST * index, 1 + REL_DIST * index]);

for i = [1, 2]
    for d = [1, 2]
        j = mod(i-1+d, 3) + 1;
        if (i < j)

            %% calculate offset of two-by-two subgame lines:
            if (i == 2) && (j == 3)
                delta = [0, -REL_DIST * index];
            elseif (i == 1) && (j == 2)
                delta = [-REL_DIST * index, REL_DIST * index * tan(pi/6)];
            else
                delta = [REL_DIST * index, REL_DIST * index * tan(pi/6)];
            end

            %% line data for offset line for two-by-two sub-game:
            line_x = zeros(length(line), 3);
            line_x(:, i) = line';
            line_x(:, j) = 1 - line';
            [x_line, y_line] = UVW_to_XY(line_x);

            %% Logic for 2-archetype subgame

            %% Plots line, equilibria points, and arrows

            %% WARNING: NOT AN EXACT SOLVER, MAY BE INCORRECT

            %% ASSUMPTION: The fitness functions are eiher entirely

            %%        equivalent or equivalent at a finite number of points
            % Determine functions
            if (i == 1 && j == 2)
                Wi = @(x) A.W1(x);
                Wj = @(x) A.W2(x);
            elseif (i == 1 && j == 3)
                Wi = @(x) A.W1(x);
                Wj = @(x) A.W3(x);
            elseif (i == 2 && j == 3)
                Wi = @(x) A.W2(x);
                Wj = @(x) A.W3(x);
            else
                print("Index error")
            end

            % Calculate difference at points
            mesh = generate_2d_grid(0.005, i, j);

            arrow_grid = cell2mat(cellfun(@(x) Wi(x)-Wj(x), mesh, 'UniformOutput', false));

            %% Casework for lines and arrows
            if (eq_throughout(arrow_grid, 10^(-6))) % Check if equivalent throughout
                if (index > 0)
                    plot(x_line+delta(1), y_line+delta(2), ':', 'LineWidth', 3, 'Color', color);
                    hold on;
                end
            else

                % Plot line - know it must be solid
                if (index > 0)
                    plot(x_line+delta(1), y_line+delta(2), '-', 'LineWidth', 3, 'Color', color);
                    hold on;
                end

                % Check if same sign throughout
                sign_val = same_sign_throughout(arrow_grid);

                if (sign_val == 1) % Goes to i
                    % i is stable
                    x = zeros(1, 3);
                    x(i) = 1;
                    [x_point, y_point] = UVW_to_XY(x);
                    plot(x_point+delta(1), y_point+delta(2), '.', 'MarkerSize', s_ms, 'Color', color);
                    hold on;

                    % j is unstable:
                    x = zeros(1, 3);
                    x(j) = 1;
                    [x_point, y_point] = UVW_to_XY(x);
                    plot(x_point+delta(1), y_point+delta(2), 'o', 'LineWidth', 1, 'MarkerSize', us_ms, 'Color', [1, 1, 1], 'MarkerEdgeColor', color, 'MarkerFaceColor', [1, 1, 1]);

                    % Plot arrow
                    plot_arrow(i, j, 1, delta, color, 1);

                elseif (sign_val == -1) % Goes to j
                    % j is stable:
                    x = zeros(1, 3);
                    x(j) = 1;
                    [x_point, y_point] = UVW_to_XY(x);
                    plot(x_point+delta(1), y_point+delta(2), '.', 'MarkerSize', s_ms, 'Color', color);
                    hold on;

                    % i is unstable:
                    x = zeros(1, 3);
                    x(i) = 1;
                    [x_point, y_point] = UVW_to_XY(x);
                    plot(x_point+delta(1), y_point+delta(2), 'o', 'LineWidth', 1, 'MarkerSize', us_ms, 'Color', [1, 1, 1], 'MarkerEdgeColor', color, 'MarkerFaceColor', [1, 1, 1]);

                    % Plot arrow
                    plot_arrow(i, j, 1, delta, color, -1);

                else % There exist internal equilibrium(s)
                    eq_lst = find_edge_eqs(mesh, arrow_grid, Wi, Wj, i, j);
                    if isempty(eq_lst)
                        disp("Didn't find any internal equilibriums.");
                    elseif length(eq_lst) > 10
                        disp("Found a suspiciously high number of internal equilibriums.");
                    else
                        arrow_lst = find_edge_arrows(eq_lst, Wi, Wj, i, j);
                        plot_arrows(i, j, eq_lst, delta, color, arrow_lst);
                        plot_eqs(i, j, eq_lst, delta, color, arrow_lst);
                    end

                end

            end
        end
    end

    dark = [0, 0, 0] + 0.5;

    %% Logic for plotting corner equilibria
    for corner = [1, 2, 3]

        x = zeros(1, 3);
        x(corner) = 1;
        [x_point, y_point] = UVW_to_XY(x);

        corner_type = corner_fp_type(A, corner);

        % only plot these if the game is not bumped out by user:
        if (effective_index == 0)
            rotate = 0;
            if isnan(corner_type)
                % neutral
                fixed_point_marker(rotate, x_point, y_point, color, dark, dark);
            else
                if (corner_type == 2)
                    % source
                    fixed_point_marker(rotate, x_point, y_point, color, [1, 1, 1], [1, 1, 1]);
                elseif (corner_type == 1)
                    % semi-source
                    fixed_point_marker(rotate, x_point, y_point, color, [1, 1, 1], dark);
                elseif (corner_type == 0)
                    % saddle
                    fixed_point_marker(rotate, x_point, y_point, color, color, [1, 1, 1]);
                elseif (corner_type == -1)
                    % semi-sink
                    fixed_point_marker(rotate, x_point, y_point, color, color, dark);
                elseif (corner_type == -2)
                    % sink
                    fixed_point_marker(rotate, x_point, y_point, color, color, color);
                end
            end
        end
    end

    marker_legend(effective_index, color);

    add_labels(labels);
end
end

function marker_legend(effective_index, color)
if (effective_index == 0)

    % order:saddle, semisource, source (these 6 could also be reversed).

    dark = [0, 0, 0] + 0.5;

    %sink
    y0 = 0.9;
    fixed_point_marker(0, -0.1, y0, color, color, color);
    text([-0.05], [y0], 'Sink', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'fontsize', 14);

    % one zero, one negative
    y0 = 0.85;
    fixed_point_marker(0, -0.1, y0, color, color, dark);
    text([-0.05], [y0], 'Semi-sink', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'fontsize', 14);

    % two zero eigenvalues
    y0 = 0.8;
    fixed_point_marker(0, -0.1, y0, color, dark, dark);
    text([-0.05], [y0], 'Neutral', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'fontsize', 14);

    % saddle
    y0 = 0.75;
    fixed_point_marker(0, -0.1, y0, color, color, [1, 1, 1]);
    text([-0.05], [y0], 'Saddle', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'fontsize', 14);


    % one zero, one positive
    y0 = 0.7;
    fixed_point_marker(0, -0.1, y0, color, [1, 1, 1], dark);
    text([-0.05], [y0], 'Semi-source', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'fontsize', 14);

    %source
    y0 = 0.65;
    fixed_point_marker(0, -0.1, y0, color, [1, 1, 1], [1, 1, 1]);
    text([-0.05], [y0], 'Source', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'fontsize', 14);

end

end

function plot_arrow(i, j, x_star, delta, color, stability)
% direction:
xDiff = zeros(1, 3);
xDiff(i) = 1;
xDiff(j) = -1;

% down-left
V = stability * xDiff(1);
U = stability * (xDiff(3) - xDiff(2)) * cos(pi/3);

arrow_length = -0.03;

if (stability < 0)
    arrow_length = 0.03;
end


len = 0.00001;

xMid = zeros(1, 3);
xMid(i) = x_star / 2 - arrow_length;
xMid(j) = 1 - x_star / 2 + arrow_length;
[xM, yM] = UVW_to_XY(xMid);

hw = 16;
hl = 20;

%% bottom up (ignore if too close to edge)
if (x_star > 0.1)
    ah = annotation('arrow', 'headStyle', 'cback1', 'HeadLength', hl, 'HeadWidth', hw);
    set(ah, 'parent', gca);
    set(ah, 'position', [xM + delta(1), yM + delta(2), len * U, len * V]);
    set(ah, 'Color', color);
end

%% top down (ignore if too close to edge)
if (x_star < 0.9)
    xMid = zeros(1, 3);
    xMid(i) = 1 - (1 - x_star) / 2 + arrow_length;
    xMid(j) = (1 - x_star) / 2 - arrow_length;

    [xM, yM] = UVW_to_XY(xMid);

    ah = annotation('arrow', 'headStyle', 'cback1', 'HeadLength', hl, 'HeadWidth', hw);
    set(ah, 'parent', gca);
    set(ah, 'position', [xM + delta(1), yM + delta(2), -len * U, -len * V]);
    set(ah, 'Color', color);
end
end

function plot_arrows(i, j, eq_lst, delta, color, arrow_lst)
pt_lst = []; % Stores the locations of all the arrows
prev = 0;
for eq = eq_lst
    pt = (eq + prev) / 2;
    pt_lst = [pt_lst, pt];
    prev = eq;
end
pt_lst = [pt_lst, (1 + prev) / 2];

for idx = 1:length(pt_lst) % Draw an arrow for each section
    if pt_lst(idx) > 0.05 && pt_lst(idx) < 0.95
        % direction:
        xDiff = zeros(1, 3);
        xDiff(i) = 1;
        xDiff(j) = -1;

        % down-left
        V = xDiff(1);
        U = (xDiff(3) - xDiff(2)) * cos(pi/3);
        len = 0.00001;
        hw = 16;
        hl = 20;
        arrow_length = -0.03;

        [xM, yM] = UVW_to_XY(insert_values(i, j, pt_lst(idx)));

        if arrow_lst(idx) > 0 % Arrow goes towards i
            ah = annotation('arrow', 'headStyle', 'cback1', 'HeadLength', hl, 'HeadWidth', hw);
            set(ah, 'parent', gca);
            set(ah, 'position', [xM + delta(1), yM + delta(2), len * U, len * V]);
            set(ah, 'Color', color);
        else % Arrow goes towards j
            ah = annotation('arrow', 'headStyle', 'cback1', 'HeadLength', hl, 'HeadWidth', hw);
            set(ah, 'parent', gca);
            set(ah, 'position', [xM + delta(1), yM + delta(2), -len * U, -len * V]);
            set(ah, 'Color', color);
        end

    end
end

end

function plot_eqs(i, j, eq_lst, delta, color, arrow_lst)
us_ms = 12;
s_ms = 50;
% Draw 0
x = zeros(1, 3);
x(i) = 1;
[x_point, y_point] = UVW_to_XY(x);
if arrow_lst(1) > 0 % stable i
    plot(x_point+delta(1), y_point+delta(2), '.', 'MarkerSize', s_ms, 'Color', color);
else % unstable i
    plot(x_point+delta(1), y_point+delta(2), 'o', 'LineWidth', 1, 'MarkerSize', us_ms, 'Color', [1, 1, 1], 'MarkerEdgeColor', color, 'MarkerFaceColor', [1, 1, 1]);
end
hold on;

% Draw 1
x = zeros(1, 3);
x(j) = 1;
[x_point, y_point] = UVW_to_XY(x);
if arrow_lst(length(arrow_lst)) < 0 % stable j
    plot(x_point+delta(1), y_point+delta(2), '.', 'MarkerSize', s_ms, 'Color', color);
else % unstable j
    plot(x_point+delta(1), y_point+delta(2), 'o', 'LineWidth', 1, 'MarkerSize', us_ms, 'Color', [1, 1, 1], 'MarkerEdgeColor', color, 'MarkerFaceColor', [1, 1, 1]);
end

% Draw all others
for idx = 1:length(eq_lst)
    [x_point, y_point] = UVW_to_XY(insert_values(i, j, eq_lst(idx)));
    if arrow_lst(idx) > 0 && arrow_lst(idx+1) < 0 % stable
        plot(x_point+delta(1), y_point+delta(2), '.', 'MarkerSize', s_ms, 'Color', color);
    else
        plot(x_point+delta(1), y_point+delta(2), 'o', 'LineWidth', 1, 'MarkerSize', us_ms, 'Color', [1, 1, 1], 'MarkerEdgeColor', color, 'MarkerFaceColor', [1, 1, 1]);
    end
end

end


function [type] = DetermineFixedPointType(x, A)

[~, lambda, ~] = hessian(x, A);

lambda1 = real(lambda(1, 1));
lambda2 = real(lambda(2, 2));

% fix numerical error for 0 eigenvalues
eps = 1e-10;

if (abs(lambda1) < eps)
    lambda1 = 0;
end

if (abs(lambda2) < eps)
    lambda2 = 0;
end


if ((lambda1 == 0) && (lambda2 == 0))
    % this is stable, but not asymptotically stable
    type = -1; % zero in the other file
elseif (((lambda1 == 0) && (lambda2 > 0)) || ((lambda1 > 0) && (lambda2 == 0)))
    % one zero, one negative
    type = 0;
elseif (((lambda1 == 0) && (lambda2 < 0)) || ((lambda1 < 0) && (lambda2 == 0)))
    % one zero, one positive == UNSTABLE
    type = 2;
elseif lambda1 * lambda2 < 0
    % this is a saddle point
    type = 3;
elseif ((lambda1 > 0) && (lambda2 > 0))
    % this is a source
    type = 4;
else
    % this is an attractor
    type = 1;
end
end

% Adapted from:
% - Salman Mashayekh (2020). Custom Marker Plot (https://www.mathworks.com/matlabcentral/fileexchange/39487-custom-marker-plot), MATLAB Central File Exchange. Retrieved July 28, 2020.

function [] = CustomMark(xData, yData, x, A, color, index)
type = DetermineFixedPointType(x, A)

color1 = color;
color2 = [0, 0, 0] + 0.3;

if (type == 4)
    % source: filled circle of white
    color1 = [1, 1, 1];
    color2 = [1, 1, 1];

elseif (type == 1)
    % attractor (sink): filled circle of color
    color1 = color;
    color2 = color;

elseif (type == 3)
    % saddle: color + white
    color1 = color;
    color2 = [1, 1, 1];

elseif (type == -1)
    % two zero eigenvalues
    color1 = [0, 0, 0] + 0.3;
    color2 = [0, 0, 0] + 0.3;

elseif (type == 0)
    % one zero eigenvalue, with stable eigenvalue
    color1 = color;
    color2 = [0, 0, 0] + 0.3;

elseif (type == 2)
    % one zero eigenvalue, with unstable eigenvalue
    color1 = color;
    color2 = [1, 1, 1];
end

if (index == 0)

    rotate = 0;

    if (x(1) == 1)
        rotate = pi / 2;
    elseif (x(2) == 1)
        rotate = -pi / 3 - pi / 2;
    elseif (x(3) == 1)
        rotate = pi / 3 + pi / 2;
    else
        if (x(1) == 0)
            rotate = pi / 2;
        elseif (x(2) == 0)
            rotate = pi / 6 + pi;
        elseif (x(3) == 0)
            rotate = -pi / 6 - pi;
        end
    end

    fixed_point_marker(rotate, xData, yData, color, color1, color2);
end
end

function [] = fixed_point_marker(rotate, xData, yData, color, color1, color2)
r = 1 / 55;
markerSize = 1;

lw = 1.5;
markerEdgeColor = color;
markerFaceColor = color;

% color-filled circle, then gray circle
for i = 1:2
    if i == 2
        r = 1 / 70; %slightly smaller
        markerFaceColor = color2;
        phi2 = rotate:0.01:(pi + rotate);
        markerDataX = r * cos(phi2);
        markerDataY = r * sin(phi2);
    else
        markerFaceColor = color1;
        phi1 = 0:0.01:(2 * pi);
        markerDataX = r * cos(phi1);
        markerDataY = r * sin(phi1);
    end

    xData = reshape(xData, length(xData), 1);
    yData = reshape(yData, length(yData), 1);
    markerDataX = markerSize * reshape(markerDataX, 1, length(markerDataX));
    markerDataY = markerSize * reshape(markerDataY, 1, length(markerDataY));

    vertX = repmat(markerDataX, length(xData), 1);
    vertX = vertX(:);
    vertY = repmat(markerDataY, length(yData), 1);
    vertY = vertY(:);

    vertX = repmat(xData, length(markerDataX), 1) + vertX;
    vertY = repmat(yData, length(markerDataY), 1) + vertY;
    faces = 0:length(xData):length(xData) * (length(markerDataY) - 1);
    faces = repmat(faces, length(xData), 1);
    faces = repmat((1:length(xData))', 1, length(markerDataY)) + faces;

    patchHndl = patch('Faces', faces, 'Vertices', [vertX, vertY]);

    if (i > 1)
        set(patchHndl, 'FaceColor', markerFaceColor, 'EdgeColor', 'none');
    else
        set(patchHndl, 'FaceColor', markerFaceColor, 'EdgeColor', markerEdgeColor, 'LineWidth', lw);
    end
end
end

function cfp_val = corner_fp_type(A, id)
%cfp_val = corner_fp_type(A,id)
%   Takes as input:
%       A - square 3-by-3 game matrix
%       id - strategy to analyze (1, 2, or 3; default: 1)
%   Returns the stability value (cfp_val) of monomorphic population of type
%   strat_num when invaded by either of other two strategies.
%   Possible values for cfp_val:
%       2 - source
%       1 - semi-source
%       0 - saddle
%      -1 - semi-sink
%      -2 - sink
%      NaN - neutral

%% WARNING: Analysis just checks if both edges points towards/away

if id == 1
    pt1 = insert_values(1, 2, 0.001);
    pt2 = insert_values(1, 3, 0.001);
    s1 = sign(A.W1(pt1)-A.W2(pt1)); % 1 if towards, -1 if away, 0 if eq
    s2 = sign(A.W1(pt2)-A.W3(pt2));
elseif id == 2
    pt1 = insert_values(2, 1, 0.001);
    pt2 = insert_values(2, 3, 0.001);
    s1 = sign(A.W2(pt1)-A.W1(pt1));
    s2 = sign(A.W2(pt2)-A.W3(pt2));
else
    pt1 = insert_values(3, 1, 0.001);
    pt2 = insert_values(3, 2, 0.001);
    s1 = sign(A.W3(pt1)-A.W1(pt1));
    s2 = sign(A.W3(pt2)-A.W2(pt2));
end

if s1 == 1 && s2 == 1 % sink
    cfp_val = -2;
elseif s1 == -1 && s2 == -1 % source
    cfp_val = 2;
elseif s1 == 0 && s2 == 0 % neutral
    cfp_val = NaN;
elseif s1 == 1 && s2 == -1 || s1 == -1 && s2 == 1 %saddle
    cfp_val = 0;
elseif s1 == 1 && s2 == 0 || s1 == 0 && s2 == 1 %semisink
    cfp_val = -1;
elseif s1 == -1 && s2 == 0 || s1 == 0 && s2 == -1 %semisource
    cfp_val = 1;
end
end

function fp_val = edge_fp_type(A, id)

%fp_val = edge_fp_type(A,id)
%   Takes as input:
%       A - square 2-by-2 game matrix
%       strat_num - strategy to analyze (1 or 2; default: 1)
%   Returns the stability value (fp_val) of monomorphic population of type
%   strat_num when invaded by the other strategy. Possible values for
%   fp_val:
%       2 - strict source
%       1 - source
%       0 - neutral
%      -1 - sink
%      -2 - strict sink

%check strat_num and number of inputs
if (nargin == 1)
    id = 1;
elseif (nargin == 0)
    fprintf('Edge analysis failed: provide a 2-by-2 game matrix A \n')
    return
else
    if (id ~= 2) && (id ~= 1)
        fprintf('Edge analysis warning: unclear strat_num, defaulting to 1')
        id = 1; %default to analyzing first strategy
    end
end


%check matrix input
sA = size(A);
if (sA(1, 1) ~= 2) || (sA(1, 2) ~= 2)
    fprintf('Edge analysis failed: Game matrix A has to be 2-by-2 \n')
    return
end

%transform matrix if analyzing strat 2
if id == 2
    A = [0, 1; 1, 0] * A * [0, 1; 1, 0]; %rotate the matrix
end

%%ANALYSIS STARTS

inv_fit = A(2, 1) - A(1, 1); %invasion fitness of type 2 against 1

if inv_fit > 0
    fp_val = 2;
elseif inv_fit < 0
    fp_val = -2;
else %neutral first-order, so source/sink not strict
    %check if type 2 can invade type 1 by drift (i.e. not strict)
    drift_fit = A(2, 2) - A(1, 2);

    if drift_fit > 0
        fp_val = 1;
    elseif drift_fit < 0
        fp_val = -1;
    else
        fp_val = 0;
    end
end
end

function vec = insert_values(i, j, g)
vec = zeros(1, 3);
vec(i) = g;
vec(j) = 1 - g;
end

function grid = generate_2d_grid(s, i, j)
ns = 0:s:1;
grid = arrayfun(@(g) insert_values(i, j, g), ns, 'UniformOutput', false);
end

function eq_val = eq_throughout(arrow_grid, tol)
eq_val = true;
for a = arrow_grid
    if abs(a) > tol
        eq_val = false;
        return;
    end
end
end

function sign_val = same_sign_throughout(arrow_grid)

first_nonzero_sign = sign(arrow_grid(find(arrow_grid ~= 0, 1, 'first')));

for a = arrow_grid
    if a ~= 0
        if sign(a) ~= first_nonzero_sign
            sign_val = 0;
            return; % Mixed signs
        end
    end
end

if first_nonzero_sign == -1
    sign_val = -1; % All negative
else
    sign_val = 1; % All positive
end
end

function eq_lst = find_edge_eqs(mesh, arrow_grid, Wi, Wj, i, j)
eq_lst = [];
f = @(x) Wi(insert_values(i, j, x)) - Wj(insert_values(i, j, x));

% Check if idx, idx + 1 contains an equilibria
for idx = 1:(length(arrow_grid) - 1)
    if (sign(arrow_grid(idx)) ~= sign(arrow_grid(idx+1))) || (abs(arrow_grid(idx)) < 1e-6) || (abs(arrow_grid(idx+1)) < 1e-6)
        try
            eq = fzero(f, [mesh{idx}(i), mesh{idx+1}(i)]);
            % Add the root to the list if it's not already included
            if ~ismember(eq, eq_lst)
                eq_lst = [eq_lst, eq];
            end
        catch
            % If fzero fails (i.e., no root in this interval), do nothing
        end
    end
end
end

function arrow_lst = find_edge_arrows(eq_lst, Wi, Wj, i, j)
arrow_lst = [];
prev = 0;
for eq = eq_lst
    val = (eq + prev) / 2;
    dir = sign(Wi(insert_values(i, j, val))-Wj(insert_values(i, j, val)));
    arrow_lst = [arrow_lst, dir];
    prev = eq;
end
val = (1 + prev) / 2;
dir = sign(Wi(insert_values(i, j, val))-Wj(insert_values(i, j, val)));
arrow_lst = [arrow_lst, dir];
end

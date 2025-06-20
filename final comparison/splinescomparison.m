clear
close all
% spline_plot_bernstein.m
% This script plots a given set of (X,Y) data using a smooth spline,
% and overlays Bézier curves (polynomial curves in the Bernstein basis)
% of various degrees based on specified control points.

% ------------------------
% User Inputs
% ------------------------
% 1) Data to curve through (cubic spline for original data)
xData = [-9.3670E-10
0.1e3-7.5640E-01
0.2e3-5.5210E+00
0.3e3-1.6630E+01
0.4e3-3.4990E+01
0.5e3-6.0510E+01
0.6e3-9.2460E+01
0.7e3-1.2980E+02
0.8e3-1.7120E+02
0.9e3-2.1540E+02
1e3-2.6100E+02
]';
yData = [3.06300E-09
-1.07800E+01
-4.08700E+01
-8.64900E+01
-1.44100E+02
-2.10800E+02
-2.84100E+02
-3.62000E+02
-4.43000E+02
-5.26000E+02
-6.09900E+02]';



% 2) Control points for Bézier curves
%    3 control points (degree 2, quadratic Bézier)
bx{1} = [0 0.502122 0.991501]*1e3;     by{1} = [0 0 -0.112514]*1e3;
%    4 control points (degree 3, cubic Bézier)
bx{2} = [0 0.334497 0.659251 0.822802]*1e3; by{2} = [0 0 -0.195121 -0.485927]*1e3;
%    5 control points (degree 4, quartic Bézier)
bx{3} = [0 0.250321 0.497298 0.607708 0.740492]*1e3; by{3} = [0 0 -0.164619 -0.386609 -0.598833]*1e3;
%    6 control points (degree 5, quintic Bézier)
bx{4} = [0 0.200155 0.39875 0.51948 0.630998 0.739303]*1e3; by{4} = [0 0 -0.115537 -0.275189 -0.441164 -0.609342]*1e3;

% Labels for plotting
labels = {"Quadratic Bézier (3 pts)", "Cubic Bézier (4 pts)", ...
          "Quartic Bézier (5 pts)", "Quintic Bézier (6 pts)"};
controlPointLabels = {"Quadratic Bézier control points", "Cubic Bézier control points", ...
                     "Quartic Bézier control points", "Quintic Bézier control points"};
polygonLabels = {"Quadratic Bézier control polygon", "Cubic Bézier control polygon", ...
                "Quartic Bézier control polygon", "Quintic Bézier control polygon"};
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], ...
          [0.9290 0.6940 0.1250], [0 1 0]};  
markers = {'^','s','o','d'};

% ------------------------
% Fine parameter vector
t = linspace(0,1,200);

% % ------------------------
% % Compute cubic spline for cleaned data
 xi = linspace(min(xData), max(xData), 200);
 yi_spline = spline(xData, yData, xi);

% ------------------------
% Prepare figure
figure; hold on; grid on;
plot(xData, yData, 'ko', 'MarkerFaceColor', 'k', 'DisplayName', 'FEM Data');
plot(xi, yi_spline, 'k-', 'LineWidth', 1.5, 'DisplayName', 'FEM Data Fitting');

% ------------------------
% Loop over Bézier sets
for k = 1:numel(bx)
    Xc = bx{k};  Yc = by{k};
    n = numel(Xc) - 1;           % Bézier degree
    Bx = zeros(size(t));
    By = zeros(size(t));
    % Bernstein basis
    for i = 0:n
        B = nchoosek(n,i) .* (1 - t).^(n - i) .* t.^i;
        Bx = Bx + B .* Xc(i+1);
        By = By + B .* Yc(i+1);
    end
    % Plot Bézier curve
    plot(Bx, By, 'Color', colors{k}, 'LineStyle', '--', ...
         'LineWidth', 1.2, 'DisplayName', labels{k});
    % Plot control points
    plot(Xc, Yc, markers{k}, 'MarkerFaceColor', colors{k}, ...
         'MarkerSize', 6, 'DisplayName', controlPointLabels{k});
    % Plot control polygon
    plot(Xc, Yc, 'Color', colors{k}, 'LineStyle', ':', ...
         'LineWidth', 1, 'DisplayName', polygonLabels{k});
end

% ------------------------
% Finalize plot
xlabel('X Coordinate (mm)'); ylabel('Y Coordinate (mm)');
title('Beam Displacement');
legend('Location','bestoutside');
hold off;

clear
close all

% 1) Data to curve through (cubic spline for original data)
% First data set from provided values
xData1 = [-9.3670E-10, 100-7.5640E-01, 200-5.5210E+00, 300-1.6630E+01, 400-3.4990E+01, 500-6.0510E+01, 600-9.2460E+01, 700-1.2980E+02, 800-1.7120E+02, 900-2.1540E+02, 1000-2.6100E+02]; % X coordinates (mm)
yData1 = [3.06300E-09, -1.07800E+01, -4.08700E+01, -8.64900E+01, -1.44100E+02, -2.10800E+02, -2.84100E+02, -3.62000E+02, -4.43000E+02, -5.26000E+02, -6.09900E+02]; % Y coordinates (mm)

% Second data set from -0.75N data (hardcoded)
xData2 = [0, 
          100 - 1.22400E+00, 
          200 - 8.74700E+00, 
          300 - 2.57900E+01, 
          400 - 5.31700E+01, 
          500 - 9.02500E+01, 
          600 - 1.35700E+02, 
          700 - 1.87800E+02, 
          800 - 2.44900E+02, 
          900 - 3.05400E+02, 
          1000 - 3.67400E+02]; % X coordinates (mm)
yData2 = [0, -1.3710E+01, -5.1260E+01, -1.0690E+02, -1.7550E+02, -2.5320E+02, -3.3700E+02, -4.2480E+02, -5.1510E+02, -6.0690E+02, -6.9950E+02]; % Y coordinates (mm)

% Third data set from -1N data (hardcoded)
xData3 = [0, 
          100 - 1.66700E+00, 
          200 - 1.16700E+01, 
          300 - 3.37500E+01, 
          400 - 6.83100E+01, 
          500 - 1.14100E+02, 
          600 - 1.69100E+02, 
          700 - 2.31300E+02, 
          800 - 2.98700E+02, 
          900 - 3.69500E+02, 
          1000 - 4.41900E+02]; % X coordinates (mm)
yData3 = [0, -1.5990E+01, -5.9030E+01, -1.2150E+02, -1.9700E+02, -2.8090E+02, -3.7020E+02, -4.6280E+02, -5.5730E+02, -6.5300E+02, -7.4910E+02]; % Y coordinates (mm)

% 2) Control points for Bézier curves (6 control points each)
% First provided Bézier curve
bx{1} = [0 0.200155 0.39875 0.51948 0.630998 0.739303]*1e3; 
by{1} = [0 0 -0.115537 -0.275189 -0.441164 -0.609342]*1e3;

% Control points for the second curve (-0.75N)
bx{2} = [0, 0.200359202850183, 0.396798388970612, 0.472192276103971, 0.554896356562057, 0.629182346388501]*1e3; 
by{2} = [0, 0, -0.149981422009211, -0.333404598046099, -0.515723300097161, -0.701466141714871]*1e3;

% Control points for the third curve (-1N)
bx{3} = [0, 0.200514159118231, 0.395120348822528, 0.440329692289411, 0.508471840277205, 0.565059930645154]*1e3;
by{3} = [0, 0, -0.170789101844841, -0.363138734402648, -0.551597881097355, -0.743492923930821]*1e3;

% Labels for plotting
dataLabels = {"Data Curve 1 (-0.5N)", "Data Curve 2 (-0.75N)", "Data Curve 3 (-1N)"};
bezierLabels = {"Bézier Curve 1 (-0.5N)", "Bézier Curve 2 (-0.75N)", "Bézier Curve 3 (-1N)"};
controlPointLabels = {"Control Points 1 (-0.5N)", "Control Points 2 (-0.75N)", "Control Points 3 (-1N)"};
polygonLabels = {"Control Polygon 1 (-0.5N)", "Control Polygon 2 (-0.75N)", "Control Polygon 3 (-1N)"};
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.9290 0.6940 0.1250]};
markers = {'o', 's', 'd'};

% ------------------------
% Fine parameter vector
t = linspace(0, 1, 200);

% ------------------------
% Prepare figure
figure; hold on; grid on;

% Plot data curves
for i = 1:3
    % Match data with the correct curve
    if i == 2
        xData = xData2;  % -0.75N data
        yData = yData2;
    elseif i == 3
        xData = xData3;  % -1N data
        yData = yData3;
    else
        xData = xData1;  % -0.5N data
        yData = yData1;
    end
    
    xi = linspace(min(xData), max(xData), 200);
    yi_spline = spline(xData, yData, xi);
    plot(xi, yi_spline, 'Color', colors{i}, 'LineWidth', 2, ...
         'DisplayName', dataLabels{i});
end

% Plot Bézier curves
for k = 1:3
    Xc = bx{k}; Yc = by{k};
    n = numel(Xc) - 1; % Bézier degree
    Bx = zeros(size(t));
    By = zeros(size(t));
    
    % Bernstein basis
    for i = 0:n
        B = nchoosek(n, i) .* (1 - t).^(n - i) .* t.^i;
        Bx = Bx + B .* Xc(i+1);
        By = By + B .* Yc(i+1);
    end
    
    % Plot Bézier curve
    plot(Bx, By, 'Color', colors{k}, 'LineStyle', '--', 'LineWidth', 1.5, ...
         'DisplayName', bezierLabels{k});
    
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
title('Comparison between FEM and spline model with different forces applied');
legend('Location', 'bestoutside');
hold off;
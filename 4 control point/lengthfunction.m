% Define control points
P0 = [0, 0];
P1 = [1/3, 0];
P2 = [2/3, 0];
P3 = [1, 0];

% Compute the length
L = bezier_cubic_length(P0, P1, P2, P3);
fprintf('The length of the cubic Bézier curve is: %.4f\n', L);

function length = bezier_cubic_length(P0, P1, P2, P3)
    % Define control point coordinates
    x0 = P0(1); y0 = P0(2);
    x1 = P1(1); y1 = P1(2);
    x2 = P2(1); y2 = P2(2);
    x3 = P3(1); y3 = P3(2);

    % Derivatives dx/dt and dy/dt as anonymous functions
    dxdt = @(t) ...
        3*(1 - t).^2 * (x1 - x0) + ...
        6*(1 - t).*t * (x2 - x1) + ...
        3*t.^2 * (x3 - x2);

    dydt = @(t) ...
        3*(1 - t).^2 * (y1 - y0) + ...
        6*(1 - t).*t * (y2 - y1) + ...
        3*t.^2 * (y3 - y2);

    % Arc length integrand
    integrand = @(t) sqrt(dxdt(t).^2 + dydt(t).^2);

    % Numerically integrate from t = 0 to t = 1
    length = integral(integrand, 0, 1);
end

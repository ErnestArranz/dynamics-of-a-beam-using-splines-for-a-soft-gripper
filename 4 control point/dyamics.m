%dynamics for a cubic spline. in this case integrals are computed
%numerically using Gaussian quadrature

clear
close all
clc
% Time span for the simulation

tspan = [0 6]; 

% Initial conditions:
% [ x1, x2, y2, x3, y3, vx1, vx2, vy2, vx3, vy3 ]
L =1;
X0 = [L/3, 2*L/3, 0, L, 0, 0, 0, 0, 0, 0 ];


% Solve the ODE system using ode45
[t, X] = ode45(@myStateSpace, tspan, X0);

% 1. Plot positions (x1, x2, y2, x3, y3)
subplot(2, 2, 1);
plot(t, X(:,1), 'b-', t, X(:,2), 'b--', t, X(:,3), 'r-', t, X(:,4), 'g-', t, X(:,5), 'g--');
xlabel('Time');
ylabel('Position');
title('Positions');
legend('x1', 'x2', 'y2', 'x3', 'y3');
grid on;

% 2. Plot velocities (vx1, vx2, vy2, vx3, vy3)
subplot(2, 2, 2);
plot(t, X(:,6), 'b-', t, X(:,7), 'b--', t, X(:,8), 'r-', t, X(:,9), 'g-', t, X(:,10), 'g--');
xlabel('Time');
ylabel('Velocity');
title('Velocities');
legend('vx1', 'vx2', 'vy2', 'vx3', 'vy3');
grid on;

% 3. Plot trajectories (x1 vs y1, x2 vs y2, x3 vs y3)
subplot(2, 2, 3);
plot(X(:,1), zeros(length(X),1), 'b-', X(:,2), X(:,3), 'r-', X(:,4), X(:,5), 'g-');
xlabel('X Position');
ylabel('Y Position');
title('Trajectories');
legend('P1', 'P2', 'P3');
grid on;
axis equal;

% 4. Plot length of the beam over time using Bézier curve length formula
subplot(2, 2, 4);
beam_length = arrayfun(@(i) bezier_length([0 0],[X(i,1), 0], [X(i,2), X(i,3)], [X(i,4), X(i,5)]), 1:length(t));
plot(t, beam_length, 'g-');
xlabel('Time');
ylabel('Beam Length');
title('Beam Length Over Time');
grid on;

% Calculate distances between control points
P0 = [0, 0]; % Origin
distances_P0P1 = zeros(length(t), 1);
distances_P1P2 = zeros(length(t), 1);
distances_P2P3 = zeros(length(t), 1);

for i = 1:length(t)
    P1 = [X(i,1), 0];
    P2 = [X(i,2), X(i,3)];
    P3 = [X(i,4), X(i,5)];
    
    % Calculate Euclidean distances
    distances_P0P1(i) = norm(P1 - P0);
    distances_P1P2(i) = norm(P2 - P1);
    distances_P2P3(i) = norm(P3 - P2);
end

% Create a new figure for distances between control points
figure;
plot(t, distances_P0P1, 'b-', t, distances_P1P2, 'r-', t, distances_P2P3, 'g-');
xlabel('Time');
ylabel('Distance');
title('Distances Between Control Points');
legend('P0-P1', 'P1-P2', 'P2-P3');
grid on;

% Calculate energies

%mechanical properties
D = 0.02;
E = 206e5;
I = (pi/64)*D^4;
rho = 1000;
A = (pi*D^2)/4; 
mu = rho*A;
L = 1;
l=L;

function length = bezier_length(P0, P1, P2, P3)
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

function dXdt = myStateSpace(t, X)
    % Extract variables from the state vector
    x_1 = X(1);
    x_2 = X(2);
    y_2 = X(3);
    x_3 = X(4);
    y_3 = X(5);
    vx1 = X(6);
    vx2 = X(7);
    vy2 = X(8);
    vx3 = X(9);
    vy3 = X(10);
    
    %mechanical properties
    D = 0.02;
    E = 206e5;
    I = (pi/64)*D^4;
    rho = 1000;
    A = (pi*D^2)/4; 
    mu = rho*A;
    L = 1;
    l=L;

    %damping coefficient
    c=3.5;

   % Define forces (if F!=0 when calculating equations of motion)
    if t >= 0 && t <= 6
        Fy3 = -1; % 
    else
        Fy3 = 0;
    end
 
ben1 = @(u) (2.*(2.*y_2 - 6.*u.*y_2 + 2.*u.*y_3 + 6.*u.^2.*y_2 - 4.*u.^2.*y_3).*(2.*x_1.*y_2 + 6.*u.^2.*x_1.*y_2 - 4.*u.^2.*x_1.*y_3 + 2.*u.^2.*x_2.*y_3 - 2.*u.^2.*x_3.*y_2 - 6.*u.*x_1.*y_2 + 2.*u.*x_1.*y_3))./(9.*(u.^2.*(2.*y_2 - 3.*u.*y_2 + u.*y_3).^2 + (x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3).^2).^3) - (8.*(3.*u.^2 - 4.*u + 1).*(x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3).*(x_1.*y_2 + 3.*u.^2.*x_1.*y_2 - 2.*u.^2.*x_1.*y_3 + u.^2.*x_2.*y_3 - u.^2.*x_3.*y_2 - 3.*u.*x_1.*y_2 + u.*x_1.*y_3).^2)./(3.*(u.^2.*(2.*y_2 - 3.*u.*y_2 + u.*y_3).^2 + (x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3).^2).^4);
 
ben2 = @(u) (4.*u.^2.*y_3.*(2.*x_1.*y_2 + 6.*u.^2.*x_1.*y_2 - 4.*u.^2.*x_1.*y_3 + 2.*u.^2.*x_2.*y_3 - 2.*u.^2.*x_3.*y_2 - 6.*u.*x_1.*y_2 + 2.*u.*x_1.*y_3))./(9.*(u.^2.*(2.*y_2 - 3.*u.*y_2 + u.*y_3).^2 + (x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3).^2).^3) + (8.*u.*(3.*u - 2).*(x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3).*(x_1.*y_2 + 3.*u.^2.*x_1.*y_2 - 2.*u.^2.*x_1.*y_3 + u.^2.*x_2.*y_3 - u.^2.*x_3.*y_2 - 3.*u.*x_1.*y_2 + u.*x_1.*y_3).^2)./(3.*(u.^2.*(2.*y_2 - 3.*u.*y_2 + u.*y_3).^2 + (x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3).^2).^4);
 
ben3 = @(u) (2.*(2.*x_1 - 6.*u.*x_1 + 6.*u.^2.*x_1 - 2.*u.^2.*x_3).*(2.*x_1.*y_2 + 6.*u.^2.*x_1.*y_2 - 4.*u.^2.*x_1.*y_3 + 2.*u.^2.*x_2.*y_3 - 2.*u.^2.*x_3.*y_2 - 6.*u.*x_1.*y_2 + 2.*u.*x_1.*y_3))./(9.*(u.^2.*(2.*y_2 - 3.*u.*y_2 + u.*y_3).^2 + (x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3).^2).^3) + (8.*u.^2.*(3.*u - 2).*(2.*y_2 - 3.*u.*y_2 + u.*y_3).*(x_1.*y_2 + 3.*u.^2.*x_1.*y_2 - 2.*u.^2.*x_1.*y_3 + u.^2.*x_2.*y_3 - u.^2.*x_3.*y_2 - 3.*u.*x_1.*y_2 + u.*x_1.*y_3).^2)./(3.*(u.^2.*(2.*y_2 - 3.*u.*y_2 + u.*y_3).^2 + (x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3).^2).^4);
 

ben4 = @(u) - (4.*u.^2.*y_2.*(2.*x_1.*y_2 + 6.*u.^2.*x_1.*y_2 - 4.*u.^2.*x_1.*y_3 + 2.*u.^2.*x_2.*y_3 - 2.*u.^2.*x_3.*y_2 - 6.*u.*x_1.*y_2 + 2.*u.*x_1.*y_3))./(9.*(u.^2.*(2.*y_2 - 3.*u.*y_2 + u.*y_3).^2 + (x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3).^2).^3) - (8.*u.^2.*(x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3).*(x_1.*y_2 + 3.*u.^2.*x_1.*y_2 - 2.*u.^2.*x_1.*y_3 + u.^2.*x_2.*y_3 - u.^2.*x_3.*y_2 - 3.*u.*x_1.*y_2 + u.*x_1.*y_3).^2)./(3.*(u.^2.*(2.*y_2 - 3.*u.*y_2 + u.*y_3).^2 + (x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3).^2).^4);
 

ben5 = @(u) (4.*u.*(x_1 - 2.*u.*x_1 + u.*x_2).*(2.*x_1.*y_2 + 6.*u.^2.*x_1.*y_2 - 4.*u.^2.*x_1.*y_3 + 2.*u.^2.*x_2.*y_3 - 2.*u.^2.*x_3.*y_2 - 6.*u.*x_1.*y_2 + 2.*u.*x_1.*y_3))./(9.*(u.^2.*(2.*y_2 - 3.*u.*y_2 + u.*y_3).^2 + (x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3).^2).^3) - (8.*u.^3.*(2.*y_2 - 3.*u.*y_2 + u.*y_3).*(x_1.*y_2 + 3.*u.^2.*x_1.*y_2 - 2.*u.^2.*x_1.*y_3 + u.^2.*x_2.*y_3 - u.^2.*x_3.*y_2 - 3.*u.*x_1.*y_2 + u.*x_1.*y_3).^2)./(3.*(u.^2.*(2.*y_2 - 3.*u.*y_2 + u.*y_3).^2 + (x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3).^2).^4);
 


str1 = @(u) (6.*(3.*(u.^2.*(2.*y_2 - 3.*u.*y_2 + u.*y_3).^2 + (x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3).^2).^(1./2) - 1).*(3.*u.^2 - 4.*u + 1).*(x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3))./(u.^2.*(2.*y_2 - 3.*u.*y_2 + u.*y_3).^2 + (x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3).^2).^(1./2);
 


str2 = @(u) -(6.*u.*(3.*u - 2).*(3.*(u.^2.*(2.*y_2 - 3.*u.*y_2 + u.*y_3).^2 + (x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3).^2).^(1./2) - 1).*(x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3))./(u.^2.*(2.*y_2 - 3.*u.*y_2 + u.*y_3).^2 + (x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3).^2).^(1./2);
 


str3 = @(u) -(6.*u.^2.*(3.*u - 2).*(3.*(u.^2.*(2.*y_2 - 3.*u.*y_2 + u.*y_3).^2 + (x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3).^2).^(1./2) - 1).*(2.*y_2 - 3.*u.*y_2 + u.*y_3))./(u.^2.*(2.*y_2 - 3.*u.*y_2 + u.*y_3).^2 + (x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3).^2).^(1./2);
 


str4 = @(u) (6.*u.^2.*(3.*(u.^2.*(2.*y_2 - 3.*u.*y_2 + u.*y_3).^2 + (x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3).^2).^(1./2) - 1).*(x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3))./(u.^2.*(2.*y_2 - 3.*u.*y_2 + u.*y_3).^2 + (x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3).^2).^(1./2);
 


str5 = @(u) (6.*u.^3.*(3.*(u.^2.*(2.*y_2 - 3.*u.*y_2 + u.*y_3).^2 + (x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3).^2).^(1./2) - 1).*(2.*y_2 - 3.*u.*y_2 + u.*y_3))./(u.^2.*(2.*y_2 - 3.*u.*y_2 + u.*y_3).^2 + (x_1 - 4.*u.*x_1 + 2.*u.*x_2 + 3.*u.^2.*x_1 - 3.*u.^2.*x_2 + u.^2.*x_3).^2).^(1./2);
 

% Using Gaussian quadrature instead of integral function
Qben1 = (-1/2)*E*I*L*gauss_quad(ben1, 0, 1);
Qben2 = (-1/2)*E*I*L*gauss_quad(ben2, 0, 1);
Qben3 = (-1/2)*E*I*L*gauss_quad(ben3, 0, 1);
Qben4 = (-1/2)*E*I*L*gauss_quad(ben4, 0, 1);
Qben5 = (-1/2)*E*I*L*gauss_quad(ben5, 0, 1);

Qstr1 = (-1/2) * E * A * L * gauss_quad(str1, 0, 1);
Qstr2 = (-1/2) * E * A * L * gauss_quad(str2, 0, 1);
Qstr3 = (-1/2) * E * A * L * gauss_quad(str3, 0, 1);
Qstr4 = (-1/2) * E * A * L * gauss_quad(str4, 0, 1);
Qstr5 = (-1/2) * E * A * L * gauss_quad(str5, 0, 1);

      % dynamics
    x_dotdot1 = (10*(10*Qben1 - 10*Qben2 + 3*Qben4 + 10*Qstr1 - 10*Qstr2 + 3*Qstr4))/(3*l*mu);

    x_dotdot2 = -(20*(5*Qben1 - 8*Qben2 + 3*Qben4 + 5*Qstr1 - 8*Qstr2 + 3*Qstr4))/(3*l*mu);

    y_dotdot2 = -(10*(Fy3 - 2*Qben3 + Qben5 - 2*Qstr3 + Qstr5))/(l*mu) -c*vy2;

    x_dotdot3 = (5*(2*Qben1 - 4*Qben2 + 3*Qben4 + 2*Qstr1 - 4*Qstr2 + 3*Qstr4))/(l*mu);

    y_dotdot3 = (2*(6*Fy3 - 5*Qben3 + 6*Qben5 - 5*Qstr3 + 6*Qstr5))/(l*mu) -c*vy3;


    % Assemble the derivative of the state vector
    dXdt = zeros(10,1);
    dXdt(1) = vx1;
    dXdt(2) = vx2;
    dXdt(3) = vy2;
    dXdt(4) = vx3;
    dXdt(5) = vy3;
    
    % The next components are the accelerations:
    dXdt(6) = x_dotdot1;
    dXdt(7) = x_dotdot2;
    dXdt(8) = y_dotdot2;
    dXdt(9) = x_dotdot3;
    dXdt(10) = y_dotdot3;
end

% Helper function for Gaussian quadrature 
function result = gauss_quad(func, a, b)
    % 7-point Gaussian quadrature for better accuracy
    % Points and weights for the standard interval [-1, 1]
    x = [-0.9491079123427585, -0.7415311855993945, -0.4058451513773972, ...
         0.0000000000000000, ...
         0.4058451513773972, 0.7415311855993945, 0.9491079123427585];
    
    w = [0.1294849661688697, 0.2797053914892766, 0.3818300505051189, ...
         0.4179591836734694, ...
         0.3818300505051189, 0.2797053914892766, 0.1294849661688697];
    
    % Transform to [a, b]
    x = 0.5*((b-a)*x + (b+a));
    w = 0.5*(b-a)*w;
    
    % Compute weighted sum
    result = 0;
    for i = 1:length(x)
        result = result + w(i) * func(x(i));
    end
end
clear
close all
clc
% Time span for the simulation

tspan = 0:0.001:1; 

% Initial conditions:
% [ x1, x2, y2, vx1, vx2, vy1, vy2 ]
% L = 1;
L = 63.7;

X0 = [L*2/5; L; 0; 0; 0.000; 0.000];
%X0 = [0.5; 0.9; 0; 0.375; 0.000; 0.000; 0.000; 0.000];  % i think you
%should consider initial bending/stretching only if the form of the spline
%without load is not a straight line

% Solve the ODE system using ode45
[t, X] = ode45(@myStateSpace, tspan, X0);

% 1. Plot positions (x1, x2, y2)
subplot(2, 2, 1);
plot(t, X(:,1), 'b-', t, X(:,2), 'b--', t, X(:,3), 'r-');
xlabel('Time');
ylabel('Position');
title('Positions');
legend('x1', 'x2','y2');
grid on;

% 2. Plot velocities (vx1, vx2, vy1, vy2)
subplot(2, 2, 2);
plot(t, X(:,4), 'b-', t, X(:,5), 'y--', t, X(:,6), 'r-');
xlabel('Time');
ylabel('Velocity');
title('Velocities');
legend('vx1', 'vx2', 'vy2');
grid on;

% 3. Plot trajectories (x1 vs y1, x2 vs y2)
subplot(2, 2, 3);
plot(X(:,1), zeros(length(X),1), 'b-', X(:,2), X(:,3), 'r-');
xlabel('X Position');
ylabel('Y Position');
title('Trajectories');
legend('P1', 'P2');
grid on;
axis equal;

% 4. Plot length of the beam over time using Bézier curve length formula
subplot(2, 2, 4);
beam_length = arrayfun(@(i) bezier_length([X(i,1), 0], [X(i,2), X(i,3)]), 1:length(t));
plot(t, beam_length, 'g-');
xlabel('Time');
ylabel('Beam Length');
title('Beam Length Over Time');
grid on;

% Calculate distances between control points
P0 = [0, 0]; % Origin (fixed point)
distances_P0P1 = zeros(length(t), 1);
distances_P1P2 = zeros(length(t), 1);

for i = 1:length(t)
    P1 = [X(i,1), 0];
    P2 = [X(i,2), X(i,3)];
    
    % Calculate Euclidean distances
    distances_P0P1(i) = norm(P1 - P0);
    distances_P1P2(i) = norm(P2 - P1);
end

% Create a new figure for distances between control points
figure();
plot(t, distances_P0P1, 'b-', t, distances_P1P2, 'r-');
xlabel('Time');
ylabel('Distance');
title('Distances Between Control Points');
legend('P0-P1', 'P1-P2');
grid on;

% Calculate energies

%mechanical properties
% D = 0.02;
% E = 206e5;
% I = (pi/64)*D^4;
% rho = 1000;
% A = (pi*D^2)/4; 
% mu = rho*A;
% L = 1;
% l=L;

EI = 4600;
EA = 1e3;
mu = 10e-6;
L = 63.7;
l=L;

%computed expression (suceptible to singularities)
%kinetic energy (no singularities)
K = (l*mu*(X(:,4)).^2)/15 + (l*mu*(X(:,5)).^2)/10 + (l*mu*(X(:,6)).^2)/10 + (l*mu*X(:,5).*X(:,4))/10;

%integral expressions
N = size(X, 1);          % number of rows
V_bending = zeros(N, 1); % preallocate
U_str = zeros(N, 1); % preallocate

for i = 1:N
int_bending = @(u) (l.*X(i,1).^2.*X(i,3).^2)./(4.*(u.^2.*X(i,3).^2 + (X(i,1) - 2.*u.*X(i,1) + u.*X(i,2)).^2).^3);
    
V_bending(i) = (1/2)*EI*L*integral(int_bending, 0, 1);
end

% for i = 1:N
%     int_str = @(u) (1-(2.*(u.^2.*X(i,3).^2 + (X(i,1) - 2.*u.*X(i,1) + u.*X(i,2)).^2).^(1/2))).^2;
%  U_str(i) = 1/2*EA*L*integral(int_str, 0, 1);
% end

% Sum of potentials and total energy
%k_curve = 1e-10;  % Penalty factor for straightness
%U_curve = 0.5 * k_curve * X(i,6)^2;
%V_total = V_bending + U_str + U_curve;  % Add to total potential
V_total = V_bending;% + U_str;
E_total = K + V_total;

% Replace the energy plot with individual figures for each energy
figure();
plot(t, K, 'b-');
xlabel('Time');
ylabel('Energy (J)');
title('Kinetic Energy');
grid on;

figure();
plot(t, V_bending, 'r-');
xlabel('Time');
ylabel('Energy (J)');
title('Bending Potential Energy');
grid on;

figure();
plot(t, U_str, 'm-');
xlabel('Time');
ylabel('Energy (J)');
title('Stretching Potential Energy');
grid on;

figure();
plot(t, V_total, 'c-');
xlabel('Time');
ylabel('Energy (J)');
title('Total Potential Energy');
grid on;

figure();
plot(t, E_total, 'g-');
xlabel('Time');
ylabel('Energy (J)');
title('Total Energy');
grid on;

% Add a new figure showing kinetic and potential energies together
figure();
plot(t, K, 'b-', t, V_total, 'r-');
xlabel('Time');
ylabel('Energy (J)');
title('Kinetic vs Potential Energy');
legend('Kinetic Energy', 'Potential Energy');
grid on;

function length = bezier_length(P1, P2)
    % Compute the length of a quadratic Bézier curve with control points:
    % P0 = (0,0), P1 = (x1,y1), P2 = (x2,y2)
    
    % Define the integrand for arc length calculation
    integrand = @(t) sqrt( ...
        (2*((P2(1) - 2*P1(1))*t + P1(1))).^2 + ...  % dx/dt component
        (2*((P2(2) - 2*P1(2))*t + P1(2))).^2 ...    % dy/dt component
    );
    
    % Numerically integrate the integrand from t=0 to t=1
    length = integral(integrand, 0, 1);
end

function dXdt = myStateSpace(t, X)
    % Extract variables from the state vector
    x_1 = X(1);
    x_2 = X(2);
    y_2 = X(3);
    vx1 = X(4);
    vx2 = X(5);
    vy2 = X(6);
    
    %mechanical properties
    % D = 0.02;
    % E = 206e5;
    % I = (pi/64)*D^4;
    % rho = 1000;
    % A = (pi*D^2)/4; 
    % mu = rho*A;
    % L = 1;
    % l=L;

EI = 4600;
EA = 1e3;
mu = 10e-6;
L = 63.7;
l=L;

    %damping coefficient
    c=30;
    %c=0;


   % Define forces (if F!=0 when calculating equations of motion)
    % if t >= 0 && t <= 3
    %     Fy2 = -0.3; % 
    %     Fx2 = -0.3;
    %     Fx1 = 0;
    % else
    %     Fy2 = 0;
    %     Fx2 = 0;
    %     Fx1 = 0;
    % end

    % Define forces (if F!=0 when calculating equations of motion)
    if t > 0 && t <= 0.5
        F_tip_magnitude = -2; % Adjust as needed
        T_x = (x_2 - x_1);     % Tangent vector components
        T_y = y_2;
        T_norm = norm([T_x T_y]);    
        % Perpendicular force components
        m = 12e-3/3; % mass per control point (equally distributed)
        g = 9.8; % gravity
        Fx2 = F_tip_magnitude * (-T_y) / T_norm;
        Fy2 = F_tip_magnitude * T_x / T_norm - m*g;
        theta = atan2(y_2, x_2 - x_1);
        Fx1 = F_tip_magnitude*sin(theta);
     
    else
        Fy2 = 0;
        Fx2 = 0;
        Fx1 = 0;
    end



    integrand1 = @(u) ...
    (x_1.*y_2.^2)./(2.*(u.^2.*y_2.^2 + (x_1 - 2.*u.*x_1 + u.*x_2).^2).^3) + (3.*x_1.^2.*y_2.^2.*(2.*u - 1).*(x_1 - 2.*u.*x_1 + u.*x_2))./(2.*(u.^2.*y_2.^2 + (x_1 - 2.*u.*x_1 + u.*x_2).^2).^4);

Qben1 = (-1/2)*EI*L*integral(integrand1, 0, 1, 'ArrayValued', true, 'AbsTol', 1e-6, 'RelTol', 1e-3);

%— integrand for Qben2
integrand2 = @(u) ...
    -(3.*u.*x_1.^2.*y_2.^2.*(x_1 - 2.*u.*x_1 + u.*x_2))./(2.*(u.^2.*y_2.^2 + (x_1 - 2.*u.*x_1 + u.*x_2).^2).^4);
Qben2 = (-1/2)*EI*L*integral(integrand2, 0, 1, 'ArrayValued', true, 'AbsTol', 1e-6, 'RelTol', 1e-3);

%— integrand for Qben3
integrand3 = @(u) ...
    (x_1.^2.*y_2)./(2.*(u.^2.*y_2.^2 + (x_1 - 2.*u.*x_1 + u.*x_2).^2).^3) - (3.*u.^2.*x_1.^2.*y_2.^3)./(2.*(u.^2.*y_2.^2 + (x_1 - 2.*u.*x_1 + u.*x_2).^2).^4);
Qben3 = (-1/2)*EI*L*integral(integrand3, 0, 1, 'ArrayValued', true, 'AbsTol', 1e-6, 'RelTol', 1e-3);

%% Strecthing Energy

    % % Define the integrand for Qstr1
    % integrand_str1 = @(u) ...
    %     -8.*(2.*u - 1).*(2.*(u.^2.*y_2.^2 + (x_1 - 2.*u.*x_1 + u.*x_2).^2).^(1./2) - 1).*(x_1 - 2.*u.*x_1 + u.*x_2);
    % Qstr1 = (-1/2) * EA * integral(integrand_str1, 0, 1, 'ArrayValued', true, 'AbsTol', 1e-6, 'RelTol', 1e-3);
    % 
    % % Define the integrand for Qstr2
    % integrand_str2 = @(u) ...
    %     8.*u.*(2.*(u.^2.*y_2.^2 + (x_1 - 2.*u.*x_1 + u.*x_2).^2).^(1./2) - 1).*(x_1 - 2.*u.*x_1 + u.*x_2);
    % Qstr2 = (-1/2) * EA * integral(integrand_str2, 0, 1, 'ArrayValued', true, 'AbsTol', 1e-6, 'RelTol', 1e-3);
    % 
    % % Define the integrand for Qstr3
    % integrand_str3 = @(u) ...
    %     8.*u.^2.*y_2.*(2.*(u.^2.*y_2.^2 + (x_1 - 2.*u.*x_1 + u.*x_2).^2).^(1./2) - 1);
    % Qstr3 = (-1/2) * EA * integral(integrand_str3, 0, 1, 'ArrayValued', true, 'AbsTol', 1e-6, 'RelTol', 1e-3);
    % 
    % 
    % 
    % % % dynamics with damping test
    % % x_dotdot1 = (6*(2*Qben1 - Qben2 + 2*Qstr1 - Qstr2))/(l*mu);
    % % x_dotdot2 = -(2*(3*Qben1 - 4*Qben2 + 3*Qstr1 - 4*Qstr2))/(l*mu); %- c*vx2
    % % y_dotdot2 = (5*(Fy2-c*vy2 + Qben3 + Qstr3))/(l*mu);
    % 
    % 
    % % force dynamics
    % x_dotdot1 = (6*(2*(Fx1- c*vx1) + 2*Qben1 - Qben2 + 2*Qstr1 - Qstr2))/(l*mu);
    % x_dotdot2 = -(2*(3*(Fx1- c*vx2) + 3*Qben1 - 4*(Qben2 + Fx2) + 3*Qstr1 - 4*Qstr2))/(l*mu); %- c*vx2
    % y_dotdot2 = (5*((Fy2- c*vy2) + Qben3 + Qstr3))/(l*mu);

   
   %% --- Length constraint enforcement ---%
    % Current Bézier derivative norm (for constraint)
    dB_du_norm = @(u) sqrt(4*( (1-u)*(x_1) + u*(x_2 - x_1) ).^2 + 4*( (1-u)*0 + u*y_2 ).^2);
    current_length = integral(dB_du_norm, 0, 1);

    % Lagrange multiplier (λ) calculation (simplified) - Penalty method
    lambda = -mu * (current_length - L) * 1e5;  % Stiff penalty term

    % Constraint forces (derivatives of the constraint w.r.t. coordinates)
    dC_dx1 = integral(@(u) 4*( (1-u).*(x_1) + u.*(x_2 - x_1) ) .* (1-u), 0, 1);
    dC_dx2 = integral(@(u) 4*( (1-u).*(x_1) + u.*(x_2 - x_1) ) .* u, 0, 1);
    dC_dy2 = integral(@(u) 4*( (1-u).*0 + u.*y_2 ) .* u, 0, 1);

    F_constraint = lambda * [dC_dx1; dC_dx2; dC_dy2];

    %--- Dynamics with constraint forces ---%
    x_dotdot1 = (6*(2*(Fx1 + Qben1) - Qben2))/(L*mu) - c*vx1 + F_constraint(1)/mu;
    x_dotdot2 = -(2*(3*(Fx1 + Qben1) - 4*(Fx2 + Qben2)))/(L*mu) - c*vx2 + F_constraint(2)/mu;
    y_dotdot2 = (5*(Fy2 + Qben3))/(L*mu) - c*vy2 + F_constraint(3)/mu;
    

    %% Assemble the derivative of the state vector
    % The first four components are the velocities:
    dXdt = zeros(6,1);
    dXdt(1) = vx1;
    dXdt(2) = vx2;
    dXdt(3) = vy2;

    
    % The next four components are the accelerations:
    dXdt(4) = x_dotdot1;
    dXdt(5) = x_dotdot2;
    dXdt(6) = y_dotdot2;
end
clear
close all
clc
% Time span for the simulation

tspan = [0 6]; 

% Initial conditions:
% [ x1, x2, y1, y2, vx1, vx2, vy1, vy2 ]
L =1;
%X0 = [L/2; L; 0; 0; 0.000; 0.000; 0.000; 0.000];
X0 = [0.5; 0.9;0; 0.375; 00.000; 0.000; 0.000; 0.00000];  % i think you
%should consider initial bending/stretching only if the form of the spline
%without load is not a straight line

% Solve the ODE system using ode45
[t, X] = ode45(@myStateSpace, tspan, X0);

% Define font sizes for figure elements (moving this up since we need it for all plots)
titleFontSize = 16;    % Size for plot titles
axisFontSize = 14;     % Size for axis numbers
labelFontSize = 14;    % Size for axis labels
legendFontSize = 12;   % Size for legends
% Let MATLAB use default line width (typically 0.5)

% Create a figure for the initial plots
figure;

% 1. Plot positions (x1, x2, y1, y2)
subplot(2, 2, 1);
plot(t, X(:,1), 'b-', t, X(:,2), 'b--', t, X(:,3), 'r-', t, X(:,4), 'r--');
xlabel('Time', 'FontSize', labelFontSize);
ylabel('Position', 'FontSize', labelFontSize);
title('Positions', 'FontSize', titleFontSize);
legend('x1', 'x2', 'y1', 'y2', 'FontSize', legendFontSize);
ax = gca;
ax.FontSize = axisFontSize;
grid on;

% 2. Plot velocities (vx1, vx2, vy1, vy2)
subplot(2, 2, 2);
plot(t, X(:,5), 'b-', t, X(:,6), 'b--', t, X(:,7), 'r-', t, X(:,8), 'r--');
xlabel('Time', 'FontSize', labelFontSize);
ylabel('Velocity', 'FontSize', labelFontSize);
title('Velocities', 'FontSize', titleFontSize);
legend('vx1', 'vx2', 'vy1', 'vy2', 'FontSize', legendFontSize);
ax = gca;
ax.FontSize = axisFontSize;
grid on;

% 3. Plot trajectories (x1 vs y1, x2 vs y2)
subplot(2, 2, 3);
plot(X(:,1), X(:,3), 'b-', X(:,2), X(:,4), 'r-');
xlabel('X Position', 'FontSize', labelFontSize);
ylabel('Y Position', 'FontSize', labelFontSize);
title('Trajectories', 'FontSize', titleFontSize);
legend('P1', 'P2', 'FontSize', legendFontSize);
ax = gca;
ax.FontSize = axisFontSize;
grid on;
axis equal;

% 4. Plot length of the beam over time using Bézier curve length formula
subplot(2, 2, 4);
beam_length = arrayfun(@(i) bezier_length([X(i,1), X(i,3)], [X(i,2), X(i,4)]), 1:length(t));
plot(t, beam_length, 'g-');
xlabel('Time', 'FontSize', labelFontSize);
ylabel('Beam Length', 'FontSize', labelFontSize);
title('Beam Length Over Time', 'FontSize', titleFontSize);
ax = gca;
ax.FontSize = axisFontSize;
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

%computed expression (suceptible to singularities)
%kinetic energy (no singularities)
K = (l*mu*X(:,5).^2)/15 + (l*mu*X(:,6).^2)/10 + (l*mu*X(:,7).^2)/15 + (l*mu*X(:,8).^2)/10 + (l*mu*X(:,5).*X(:,6))/10 + (l*mu*X(:,7).*X(:,8))/10;

%Uben
%V_bending = (E.*I.*((3.*l.*(atan((2.*X(:,1).^2 - 3.*X(:,1).*X(:,2) + X(:,2).^2 + 2.*X(:,3).^2 - 3.*X(:,3).*X(:,4) + X(:,4).^2)./(X(:,1).*X(:,4) - X(:,2).*X(:,3))) - atan((- 2.*X(:,1).^2 + X(:,2).*X(:,1) - 2.*X(:,3).^2 + X(:,4).*X(:,3))./(X(:,1).*X(:,4) - X(:,2).*X(:,3)))).*(4.*X(:,1).^2 - 4.*X(:,1).*X(:,2) + X(:,2).^2 + 4.*X(:,3).^2 - 4.*X(:,3).*X(:,4) + X(:,4).^2).^2)./(32.*(X(:,1).*X(:,4) - X(:,2).*X(:,3)).^3) + (l.*(24.*X(:,1).^6 - 108.*X(:,1).^5.*X(:,2) + 198.*X(:,1).^4.*X(:,2).^2 + 72.*X(:,1).^4.*X(:,3).^2 - 108.*X(:,1).^4.*X(:,3).*X(:,4) + 46.*X(:,1).^4.*X(:,4).^2 - 189.*X(:,1).^3.*X(:,2).^3 - 216.*X(:,1).^3.*X(:,2).*X(:,3).^2 + 304.*X(:,1).^3.*X(:,2).*X(:,3).*X(:,4) - 123.*X(:,1).^3.*X(:,2).*X(:,4).^2 + 99.*X(:,1).^2.*X(:,2).^4 + 244.*X(:,1).^2.*X(:,2).^2.*X(:,3).^2 - 321.*X(:,1).^2.*X(:,2).^2.*X(:,3).*X(:,4) + 122.*X(:,1).^2.*X(:,2).^2.*X(:,4).^2 + 72.*X(:,1).^2.*X(:,3).^4 - 216.*X(:,1).^2.*X(:,3).^3.*X(:,4) + 244.*X(:,1).^2.*X(:,3).^2.*X(:,4).^2 - 123.*X(:,1).^2.*X(:,3).*X(:,4).^3 + 23.*X(:,1).^2.*X(:,4).^4 - 27.*X(:,1).*X(:,2).^5 - 123.*X(:,1).*X(:,2).^3.*X(:,3).^2 + 152.*X(:,1).*X(:,2).^3.*X(:,3).*X(:,4) - 54.*X(:,1).*X(:,2).^3.*X(:,4).^2 - 108.*X(:,1).*X(:,2).*X(:,3).^4 + 304.*X(:,1).*X(:,2).*X(:,3).^3.*X(:,4) - 321.*X(:,1).*X(:,2).*X(:,3).^2.*X(:,4).^2 + 152.*X(:,1).*X(:,2).*X(:,3).*X(:,4).^3 - 27.*X(:,1).*X(:,2).*X(:,4).^4 + 3.*X(:,2).^6 + 23.*X(:,2).^4.*X(:,3).^2 - 27.*X(:,2).^4.*X(:,3).*X(:,4) + 9.*X(:,2).^4.*X(:,4).^2 + 46.*X(:,2).^2.*X(:,3).^4 - 123.*X(:,2).^2.*X(:,3).^3.*X(:,4) + 122.*X(:,2).^2.*X(:,3).^2.*X(:,4).^2 - 54.*X(:,2).^2.*X(:,3).*X(:,4).^3 + 9.*X(:,2).^2.*X(:,4).^4 + 24.*X(:,3).^6 - 108.*X(:,3).^5.*X(:,4) + 198.*X(:,3).^4.*X(:,4).^2 - 189.*X(:,3).^3.*X(:,4).^3 + 99.*X(:,3).^2.*X(:,4).^4 - 27.*X(:,3).*X(:,4).^5 + 3.*X(:,4).^6))./(32.*(X(:,1).*X(:,4) - X(:,2).*X(:,3)).^2.*(X(:,1).^2 - 2.*X(:,1).*X(:,2) + X(:,2).^2 + X(:,3).^2 - 2.*X(:,3).*X(:,4) + X(:,4).^2).^2) + (l.*(24.*X(:,1).^6 - 36.*X(:,1).^5.*X(:,2) + 18.*X(:,1).^4.*X(:,2).^2 + 72.*X(:,1).^4.*X(:,3).^2 - 36.*X(:,1).^4.*X(:,3).*X(:,4) + 10.*X(:,1).^4.*X(:,4).^2 - 3.*X(:,1).^3.*X(:,2).^3 - 72.*X(:,1).^3.*X(:,2).*X(:,3).^2 + 16.*X(:,1).^3.*X(:,2).*X(:,3).*X(:,4) - 5.*X(:,1).^3.*X(:,2).*X(:,4).^2 + 28.*X(:,1).^2.*X(:,2).^2.*X(:,3).^2 + X(:,1).^2.*X(:,2).^2.*X(:,3).*X(:,4) + 72.*X(:,1).^2.*X(:,3).^4 - 72.*X(:,1).^2.*X(:,3).^3.*X(:,4) + 28.*X(:,1).^2.*X(:,3).^2.*X(:,4).^2 - 5.*X(:,1).^2.*X(:,3).*X(:,4).^3 - 5.*X(:,1).*X(:,2).^3.*X(:,3).^2 - 36.*X(:,1).*X(:,2).*X(:,3).^4 + 16.*X(:,1).*X(:,2).*X(:,3).^3.*X(:,4) + X(:,1).*X(:,2).*X(:,3).^2.*X(:,4).^2 + 10.*X(:,2).^2.*X(:,3).^4 - 5.*X(:,2).^2.*X(:,3).^3.*X(:,4) + 24.*X(:,3).^6 - 36.*X(:,3).^5.*X(:,4) + 18.*X(:,3).^4.*X(:,4).^2 - 3.*X(:,3).^3.*X(:,4).^3))./(32.*(X(:,1).*X(:,4) - X(:,2).*X(:,3)).^2.*(X(:,1).^2 + X(:,3).^2).^2)))./2;
 
%Ustretching
%U_str = -(A.*E.*((8.*X(:,1).*X(:,2))./3 + (8.*X(:,3).*X(:,4))./3 - (8.*X(:,1).^2)./3 - (8.*X(:,2).^2)./3 - (8.*X(:,3).^2)./3 - (8.*X(:,4).^2)./3 + ((X(:,1).^2 + X(:,3).^2).^(1/2).*(16.*X(:,1).^6 - 24.*X(:,1).^5.*X(:,2) + 12.*X(:,1).^4.*X(:,2).^2 + 48.*X(:,1).^4.*X(:,3).^2 - 24.*X(:,1).^4.*X(:,3).*X(:,4) + 10.*X(:,1).^4.*X(:,4).^2 + 8.*X(:,1).^4 - 2.*X(:,1).^3.*X(:,2).^3 - 48.*X(:,1).^3.*X(:,2).*X(:,3).^2 + 4.*X(:,1).^3.*X(:,2).*X(:,3).*X(:,4) - 5.*X(:,1).^3.*X(:,2).*X(:,4).^2 - 12.*X(:,1).^3.*X(:,2) + 22.*X(:,1).^2.*X(:,2).^2.*X(:,3).^2 + 4.*X(:,1).^2.*X(:,2).^2.*X(:,3).*X(:,4) + 6.*X(:,1).^2.*X(:,2).^2 + 48.*X(:,1).^2.*X(:,3).^4 - 48.*X(:,1).^2.*X(:,3).^3.*X(:,4) + 22.*X(:,1).^2.*X(:,3).^2.*X(:,4).^2 + 16.*X(:,1).^2.*X(:,3).^2 - 5.*X(:,1).^2.*X(:,3).*X(:,4).^3 - 12.*X(:,1).^2.*X(:,3).*X(:,4) + 2.*X(:,1).^2.*X(:,4).^2 - 5.*X(:,1).*X(:,2).^3.*X(:,3).^2 - X(:,1).*X(:,2).^3 - 24.*X(:,1).*X(:,2).*X(:,3).^4 + 4.*X(:,1).*X(:,2).*X(:,3).^3.*X(:,4) + 4.*X(:,1).*X(:,2).*X(:,3).^2.*X(:,4).^2 - 12.*X(:,1).*X(:,2).*X(:,3).^2 + 8.*X(:,1).*X(:,2).*X(:,3).*X(:,4) - X(:,1).*X(:,2).*X(:,4).^2 + 10.*X(:,2).^2.*X(:,3).^4 - 5.*X(:,2).^2.*X(:,3).^3.*X(:,4) + 2.*X(:,2).^2.*X(:,3).^2 - X(:,2).^2.*X(:,3).*X(:,4) + 16.*X(:,3).^6 - 24.*X(:,3).^5.*X(:,4) + 12.*X(:,3).^4.*X(:,4).^2 + 8.*X(:,3).^4 - 2.*X(:,3).^3.*X(:,4).^3 - 12.*X(:,3).^3.*X(:,4) + 6.*X(:,3).^2.*X(:,4).^2 - X(:,3).*X(:,4).^3))./(4.*X(:,1).^2 - 4.*X(:,1).*X(:,2) + X(:,2).^2 + 4.*X(:,3).^2 - 4.*X(:,3).*X(:,4) + X(:,4).^2).^2 + ((X(:,1).^2 - 2.*X(:,1).*X(:,2) + X(:,2).^2 + X(:,3).^2 - 2.*X(:,3).*X(:,4) + X(:,4).^2).^(1/2).*(16.*X(:,1).^6 - 72.*X(:,1).^5.*X(:,2) + 132.*X(:,1).^4.*X(:,2).^2 + 48.*X(:,1).^4.*X(:,3).^2 - 72.*X(:,1).^4.*X(:,3).*X(:,4) + 34.*X(:,1).^4.*X(:,4).^2 + 8.*X(:,1).^4 - 126.*X(:,1).^3.*X(:,2).^3 - 144.*X(:,1).^3.*X(:,2).*X(:,3).^2 + 196.*X(:,1).^3.*X(:,2).*X(:,3).*X(:,4) - 87.*X(:,1).^3.*X(:,2).*X(:,4).^2 - 20.*X(:,1).^3.*X(:,2) + 66.*X(:,1).^2.*X(:,2).^4 + 166.*X(:,1).^2.*X(:,2).^2.*X(:,3).^2 - 204.*X(:,1).^2.*X(:,2).^2.*X(:,3).*X(:,4) + 83.*X(:,1).^2.*X(:,2).^2.*X(:,4).^2 + 18.*X(:,1).^2.*X(:,2).^2 + 48.*X(:,1).^2.*X(:,3).^4 - 144.*X(:,1).^2.*X(:,3).^3.*X(:,4) + 166.*X(:,1).^2.*X(:,3).^2.*X(:,4).^2 + 16.*X(:,1).^2.*X(:,3).^2 - 87.*X(:,1).^2.*X(:,3).*X(:,4).^3 - 20.*X(:,1).^2.*X(:,3).*X(:,4) + 17.*X(:,1).^2.*X(:,4).^4 + 6.*X(:,1).^2.*X(:,4).^2 - 18.*X(:,1).*X(:,2).^5 - 87.*X(:,1).*X(:,2).^3.*X(:,3).^2 + 98.*X(:,1).*X(:,2).^3.*X(:,3).*X(:,4) - 36.*X(:,1).*X(:,2).^3.*X(:,4).^2 - 7.*X(:,1).*X(:,2).^3 - 72.*X(:,1).*X(:,2).*X(:,3).^4 + 196.*X(:,1).*X(:,2).*X(:,3).^3.*X(:,4) - 204.*X(:,1).*X(:,2).*X(:,3).^2.*X(:,4).^2 - 20.*X(:,1).*X(:,2).*X(:,3).^2 + 98.*X(:,1).*X(:,2).*X(:,3).*X(:,4).^3 + 24.*X(:,1).*X(:,2).*X(:,3).*X(:,4) - 18.*X(:,1).*X(:,2).*X(:,4).^4 - 7.*X(:,1).*X(:,2).*X(:,4).^2 + 2.*X(:,2).^6 + 17.*X(:,2).^4.*X(:,3).^2 - 18.*X(:,2).^4.*X(:,3).*X(:,4) + 6.*X(:,2).^4.*X(:,4).^2 + X(:,2).^4 + 34.*X(:,2).^2.*X(:,3).^4 - 87.*X(:,2).^2.*X(:,3).^3.*X(:,4) + 83.*X(:,2).^2.*X(:,3).^2.*X(:,4).^2 + 6.*X(:,2).^2.*X(:,3).^2 - 36.*X(:,2).^2.*X(:,3).*X(:,4).^3 - 7.*X(:,2).^2.*X(:,3).*X(:,4) + 6.*X(:,2).^2.*X(:,4).^4 + 2.*X(:,2).^2.*X(:,4).^2 + 16.*X(:,3).^6 - 72.*X(:,3).^5.*X(:,4) + 132.*X(:,3).^4.*X(:,4).^2 + 8.*X(:,3).^4 - 126.*X(:,3).^3.*X(:,4).^3 - 20.*X(:,3).^3.*X(:,4) + 66.*X(:,3).^2.*X(:,4).^4 + 18.*X(:,3).^2.*X(:,4).^2 - 18.*X(:,3).*X(:,4).^5 - 7.*X(:,3).*X(:,4).^3 + 2.*X(:,4).^6 + X(:,4).^4))./(4.*X(:,1).^2 - 4.*X(:,1).*X(:,2) + X(:,2).^2 + 4.*X(:,3).^2 - 4.*X(:,3).*X(:,4) + X(:,4).^2).^2 - (log((- 2.*X(:,1).^2 + X(:,2).*X(:,1) - 2.*X(:,3).^2 + X(:,4).*X(:,3))./((2.*X(:,1) - X(:,2)).^2 + (2.*X(:,3) - X(:,4)).^2).^(1/2) + (X(:,1).^2 + X(:,3).^2).^(1/2)).*(X(:,1).*X(:,4) - X(:,2).*X(:,3)).^2.*(3.*X(:,1).^2.*X(:,4).^2 + 4.*X(:,1).^2 - 6.*X(:,1).*X(:,2).*X(:,3).*X(:,4) - 4.*X(:,1).*X(:,2) + 3.*X(:,2).^2.*X(:,3).^2 + X(:,2).^2 + 4.*X(:,3).^2 - 4.*X(:,3).*X(:,4) + X(:,4).^2))./(((2.*X(:,1) - X(:,2)).^2 + (2.*X(:,3) - X(:,4)).^2).^(1/2).*(4.*X(:,1).^2 - 4.*X(:,1).*X(:,2) + X(:,2).^2 + 4.*X(:,3).^2 - 4.*X(:,3).*X(:,4) + X(:,4).^2).^2) + (log((2.*X(:,1).^2 - 3.*X(:,1).*X(:,2) + X(:,2).^2 + 2.*X(:,3).^2 - 3.*X(:,3).*X(:,4) + X(:,4).^2)./((2.*X(:,1) - X(:,2)).^2 + (2.*X(:,3) - X(:,4)).^2).^(1/2) + (X(:,1).^2 - 2.*X(:,1).*X(:,2) + X(:,2).^2 + X(:,3).^2 - 2.*X(:,3).*X(:,4) + X(:,4).^2).^(1/2)).*(X(:,1).*X(:,4) - X(:,2).*X(:,3)).^2.*(3.*X(:,1).^2.*X(:,4).^2 + 4.*X(:,1).^2 - 6.*X(:,1).*X(:,2).*X(:,3).*X(:,4) - 4.*X(:,1).*X(:,2) + 3.*X(:,2).^2.*X(:,3).^2 + X(:,2).^2 + 4.*X(:,3).^2 - 4.*X(:,3).*X(:,4) + X(:,4).^2))./(((2.*X(:,1) - X(:,2)).^2 + (2.*X(:,3) - X(:,4)).^2).^(1/2).*(4.*X(:,1).^2 - 4.*X(:,1).*X(:,2) + X(:,2).^2 + 4.*X(:,3).^2 - 4.*X(:,3).*X(:,4) + X(:,4).^2).^2)))./2;


%integral expressions
N = size(X, 1);          % number of rows
V_bending = zeros(N, 1); % preallocate
U_str = zeros(N, 1); % preallocate

for i = 1:N
int_bending = @(u) ((X(i,1).*X(i,4) - X(i,2).*X(i,3)).^2)./(4.*((X(i,1) - 2.*u.*X(i,1) + u.*X(i,2)).^2 + (X(i,3) - 2.*u.*X(i,3) + u.*X(i,4)).^2).^3);
    
V_bending(i) = (1/2)*E*I*L*integral(int_bending, 0, 1);
end

for i = 1:N
    int_str = @(u) 2.*((X(i,1) - 2.*u.*X(i,1) + u.*X(i,2)).^2 + (X(i,3) - 2.*u.*X(i,3) + u.*X(i,4)).^2).^(1/2).*(2.*((X(i,1) - 2.*u.*X(i,1) + u.*X(i,2)).^2 + (X(i,3) - 2.*u.*X(i,3) + u.*X(i,4)).^2).^(1/2) - 1).^2;
 U_str(i) = 1/2*E*A*integral(int_str, 0, 1);
end


% Sum of potentials and total energy
V_total = V_bending + U_str;
E_total = K + V_total;

% Replace the energy plot with individual figures for each energy
figure();
plot(t, K, 'b-');
xlabel('Time', 'FontSize', labelFontSize);
ylabel('Energy (J)', 'FontSize', labelFontSize);
title('Kinetic Energy', 'FontSize', titleFontSize);
ax = gca;
ax.FontSize = axisFontSize;
grid on;

figure();
plot(t, V_bending, 'r-');
xlabel('Time', 'FontSize', labelFontSize);
ylabel('Energy (J)', 'FontSize', labelFontSize);
title('Bending Potential Energy', 'FontSize', titleFontSize);
ax = gca;
ax.FontSize = axisFontSize;
grid on;

figure();
plot(t, U_str, 'm-');
xlabel('Time', 'FontSize', labelFontSize);
ylabel('Energy (J)', 'FontSize', labelFontSize);
title('Stretching Potential Energy', 'FontSize', titleFontSize);
ax = gca;
ax.FontSize = axisFontSize;
grid on;

figure();
plot(t, V_total, 'c-');
xlabel('Time', 'FontSize', labelFontSize);
ylabel('Energy (J)', 'FontSize', labelFontSize);
title('Total Potential Energy', 'FontSize', titleFontSize);
ax = gca;
ax.FontSize = axisFontSize;
grid on;

figure();
plot(t, E_total, 'g-');
xlabel('Time', 'FontSize', labelFontSize);
ylabel('Energy (J)', 'FontSize', labelFontSize);
title('Total Energy', 'FontSize', titleFontSize);
ax = gca;
ax.FontSize = axisFontSize;
grid on;

% Add a new figure showing kinetic and potential energies together
figure();
plot(t, K, 'b-', t, V_total, 'r-');
xlabel('Time', 'FontSize', labelFontSize);
ylabel('Energy (J)', 'FontSize', labelFontSize);
title('Kinetic vs Potential Energy', 'FontSize', titleFontSize);
legend('Kinetic Energy', 'Potential Energy', 'FontSize', legendFontSize);
ax = gca;
ax.FontSize = axisFontSize;
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
  % y_1 = X(3);
   y_1 = 0; %fixed end, cantilever
    y_2 = X(4);
    vx1 = X(5);
    vx2 = X(6);
   % vy1 = X(7);
    vy1 = 0; %fixed end, cantilever
    vy2 = X(8);
    
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
    c=2;
    %c=0;


   % Define forces (if F!=0 when calculating equations of motion)
    if t >= 0 && t <= 3
        Fy2 = -0.0; % 
    else
        Fy2 = 0;
    end

    integrand1 = @(u) ...
    (y_2 .* (4 .* x_1 .* y_2 - 4 .* x_2 .* y_1)) ...
      ./ (8 .* ((x_1 - 2 .* u .* x_1 + u .* x_2).^2 + (y_1 - 2 .* u .* y_1 + u .* y_2).^2).^3) ...
    + (3 .* (x_1 .* y_2 - x_2 .* y_1).^2 .* (2 .* u - 1) .* (x_1 - 2 .* u .* x_1 + u .* x_2)) ...
      ./ (2 .* ((x_1 - 2 .* u .* x_1 + u .* x_2).^2 + (y_1 - 2 .* u .* y_1 + u .* y_2).^2).^4);

Qben1 = (-1/2)*E*I*L*integral(integrand1, 0, 1, 'ArrayValued', true);

%— integrand for Qben2
integrand2 = @(u) ...
    -(y_1 .* (4 .* x_1 .* y_2 - 4 .* x_2 .* y_1)) ...
      ./ (8 .* ((x_1 - 2 .* u .* x_1 + u .* x_2).^2 + (y_1 - 2 .* u .* y_1 + u .* y_2).^2).^3) ...
    - (3 .* u .* (x_1 .* y_2 - x_2 .* y_1).^2 .* (x_1 - 2 .* u .* x_1 + u .* x_2)) ...
      ./ (2 .* ((x_1 - 2 .* u .* x_1 + u .* x_2).^2 + (y_1 - 2 .* u .* y_1 + u .* y_2).^2).^4);
Qben2 = (-1/2)*E*I*L*integral(integrand2, 0, 1, 'ArrayValued', true);

%— integrand for Qben3
integrand3 = @(u) ...
    (3 .* (x_1 .* y_2 - x_2 .* y_1).^2 .* (2 .* u - 1) .* (y_1 - 2 .* u .* y_1 + u .* y_2)) ...
      ./ (2 .* ((x_1 - 2 .* u .* x_1 + u .* x_2).^2 + (y_1 - 2 .* u .* y_1 + u .* y_2).^2).^4) ...
    - (x_2 .* (4 .* x_1 .* y_2 - 4 .* x_2 .* y_1)) ...
      ./ (8 .* ((x_1 - 2 .* u .* x_1 + u .* x_2).^2 + (y_1 - 2 .* u .* y_1 + u .* y_2).^2).^3);
Qben3 = (-1/2)*E*I*L*integral(integrand3, 0, 1, 'ArrayValued', true);

%— integrand for Qben4
integrand4 = @(u) ...
    (x_1 .* (4 .* x_1 .* y_2 - 4 .* x_2 .* y_1)) ...
      ./ (8 .* ((x_1 - 2 .* u .* x_1 + u .* x_2).^2 + (y_1 - 2 .* u .* y_1 + u .* y_2).^2).^3) ...
    - (3 .* u .* (x_1 .* y_2 - x_2 .* y_1).^2 .* (y_1 - 2 .* u .* y_1 + u .* y_2)) ...
      ./ (2 .* ((x_1 - 2 .* u .* x_1 + u .* x_2).^2 + (y_1 - 2 .* u .* y_1 + u .* y_2).^2).^4);
Qben4 = (-1/2)*E*I*L*integral(integrand4, 0, 1, 'ArrayValued', true);

    % Define the integrand for Qstr1
    % integrand_str1 = @(u) ...
    %     (-8 .* (2 * sqrt((x_1 - 2 .* u .* x_1 + u .* x_2).^2 + (y_1 - 2 .* u .* y_1 + u .* y_2).^2) - 1) ...
    %     .* (2 .* u - 1) .* (x_1 - 2 .* u .* x_1 + u .* x_2));
    % Qstr1 = (-1/2) * E * A * integral(integrand_str1, 0, 1, 'ArrayValued', true);
    % 
    % % Define the integrand for Qstr2
    % integrand_str2 = @(u) ...
    %     (8 .* u .* (2 * sqrt((x_1 - 2 .* u .* x_1 + u .* x_2).^2 + (y_1 - 2 .* u .* y_1 + u .* y_2).^2) - 1) ...
    %     .* (x_1 - 2 .* u .* x_1 + u .* x_2));
    % Qstr2 = (-1/2) * E * A * integral(integrand_str2, 0, 1, 'ArrayValued', true);
    % 
    % % Define the integrand for Qstr3
    % integrand_str3 = @(u) ...
    %     (-8 .* (2 * sqrt((x_1 - 2 .* u .* x_1 + u .* x_2).^2 + (y_1 - 2 .* u .* y_1 + u .* y_2).^2) - 1) ...
    %     .* (2 .* u - 1) .* (y_1 - 2 .* u .* y_1 + u .* y_2));
    % Qstr3 = (-1/2) * E * A * integral(integrand_str3, 0, 1, 'ArrayValued', true);
    % 
    % % Define the integrand for Qstr4
    % integrand_str4 = @(u) ...
    %     (8 .* u .* (2 * sqrt((x_1 - 2 .* u .* x_1 + u .* x_2).^2 + (y_1 - 2 .* u .* y_1 + u .* y_2).^2) - 1) ...
    %     .* (y_1 - 2 .* u .* y_1 + u .* y_2));
    % Qstr4 = (-1/2) * E * A * integral(integrand_str4, 0, 1, 'ArrayValued', true);

%%==========stretching simplified
integrand_str1simp =  @(u) -(4.*(2.*((x_1 - 2.*u.*x_1 + u.*x_2).^2 + (y_1 - 2.*u.*y_1 + u.*y_2).^2).^(1./2) - 1).*(2.*u - 1).*(x_1 - 2.*u.*x_1 + u.*x_2))./((x_1 - 2.*u.*x_1 + u.*x_2).^2 + (y_1 - 2.*u.*y_1 + u.*y_2).^2).^(1./2);


integrand_str2simp =  @(u) (4.*u.*(2.*((x_1 - 2.*u.*x_1 + u.*x_2).^2 + (y_1 - 2.*u.*y_1 + u.*y_2).^2).^(1./2) - 1).*(x_1 - 2.*u.*x_1 + u.*x_2))./((x_1 - 2.*u.*x_1 + u.*x_2).^2 + (y_1 - 2.*u.*y_1 + u.*y_2).^2).^(1./2);


integrand_str3simp =  @(u) -(4.*(2.*((x_1 - 2.*u.*x_1 + u.*x_2).^2 + (y_1 - 2.*u.*y_1 + u.*y_2).^2).^(1./2) - 1).*(2.*u - 1).*(y_1 - 2.*u.*y_1 + u.*y_2))./((x_1 - 2.*u.*x_1 + u.*x_2).^2 + (y_1 - 2.*u.*y_1 + u.*y_2).^2).^(1./2);


integrand_str4simp =  @(u) (4.*u.*(2.*((x_1 - 2.*u.*x_1 + u.*x_2).^2 + (y_1 - 2.*u.*y_1 + u.*y_2).^2).^(1./2) - 1).*(y_1 - 2.*u.*y_1 + u.*y_2))./((x_1 - 2.*u.*x_1 + u.*x_2).^2 + (y_1 - 2.*u.*y_1 + u.*y_2).^2).^(1./2);

 Qstr1simp = (-1/2) * E * A * L*integral(integrand_str1simp, 0, 1, 'ArrayValued', true);
 Qstr2simp = (-1/2) * E * A * L*integral(integrand_str2simp, 0, 1, 'ArrayValued', true);
 Qstr3simp= (-1/2) * E * A * L*integral(integrand_str3simp, 0, 1, 'ArrayValued', true);
 Qstr4simp = (-1/2) * E * A *L* integral(integrand_str4simp, 0, 1, 'ArrayValued', true);




    % dynamics no force
    % x_dotdot1 = (60*(2*Qben1 - Qben2 + 2*Qstr1 - Qstr2))/pi; 
    % x_dotdot2 = -(20*(3*Qben1 - 4*Qben2 + 3*Qstr1 - 4*Qstr2))/pi; 
    % y_dotdot1 = (60*(2*Qben3 - Qben4 + 2*Qstr3 - Qstr4))/pi; 
    % y_dotdot2 = -(20*(3*Qben3 - 4*Qben4 + 3*Qstr3 - 4*Qstr4))/pi;
    %generic expressions
    % x_dotdot1 =           (6*(2*Qben1 - Qben2 + 2*Qstr1 - Qstr2))/(l*mu);
    % x_dotdot2 =    -(2*(3*Qben1 - 4*Qben2 + 3*Qstr1 - 4*Qstr2))/(l*mu);
    % y_dotdot1 =  -(6*(Fy2 - 2*Qben3 + Qben4 - 2*Qstr3 + Qstr4))/(l*mu);
    % y_dotdot2 =(2*(4*Fy2 - 3*Qben3 + 4*Qben4 - 3*Qstr3 + 4*Qstr4))/(l*mu);


    %  dynamics with damping test
    x_dotdot1 = (6*(2*Qben1 - Qben2 + 2*Qstr1simp - Qstr2simp))/(l*mu) - c*vx1;
    x_dotdot2 = -(2*(3*Qben1 - 4*Qben2 + 3*Qstr1simp - 4*Qstr2simp))/(l*mu) - c*vx2;
    y_dotdot1 =0;%cantilever %-(6*(Fy2 - 2*Qben3 + Qben4 - 2*Qstr3 + Qstr4))/(l*mu) - c*vy1;
    y_dotdot2 = (2*(4*Fy2 - 3*Qben3 + 4*Qben4 - 3*Qstr3simp + 4*Qstr4simp))/(l*mu) - c*vy2;

      % dynamics with damping test
    % x_dotdot1 = (6*(2*Qben1 - Qben2 + 2*Qstr1 - Qstr2))/(l*mu); %- c*vx1;
    % x_dotdot2 = -(2*(3*Qben1 - 4*Qben2 + 3*Qstr1 - 4*Qstr2))/(l*mu); %- c*vx2;
    % y_dotdot1 =-(6*(Fy2 - 2*Qben3 + Qben4 - 2*Qstr3 + Qstr4))/(l*mu); %- c*vy1;
    % y_dotdot2 = (2*(4*(Fy2- c*vy2) - 3*Qben3 + 4*Qben4 - 3*Qstr3 + 4*Qstr4))/(l*mu) ;

%cantilever dynamics 3 terms deleted
    % x_dotdot1 = (6*(2*Qben1 - Qben2 + 2*Qstr1 - Qstr2))/(l*mu); %- c*vx1;
    % x_dotdot2 = -(2*(3*Qben1 - 4*Qben2 + 3*Qstr1 - 4*Qstr2))/(l*mu); %- c*vx2;
    % y_dotdot1 =0;
    % y_dotdot2 = (2*(4*(Fy2- c*vy2) - 3*0 + 4*Qben4 - 3*0 + 4*Qstr4))/(l*mu) ;
    % % 


    % Assemble the derivative of the state vector
    % The first four components are the velocities:
    dXdt = zeros(8,1);
    dXdt(1) = vx1;
    dXdt(2) = vx2;
    dXdt(3) = vy1;
    dXdt(4) = vy2;
    
    % The next four components are the accelerations:
    dXdt(5) = x_dotdot1;
    dXdt(6) = x_dotdot2;
    dXdt(7) = y_dotdot1;
    dXdt(8) = y_dotdot2;
end
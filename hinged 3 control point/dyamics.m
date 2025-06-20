
clear
close
% Time span for the simulation



tspan = [0 10]; % from t = 0 to t = 10

% Initial conditions:
% [ x1, x2, y1, y2, vx1, vx2, vy1, vy2 ]
L =1;
X0 = [L/2; L-1e-6; -1e-6; -2e-6; 0.000; 0.000; 0.000; 0.000];  % example initial conditions

% Solve the ODE system using ode45
[t, X] = ode45(@myStateSpace, tspan, X0);

% 1. Plot positions (x1, x2, y1, y2)
subplot(2, 2, 1);
plot(t, X(:,1), 'b-', t, X(:,2), 'b--', t, X(:,3), 'r-', t, X(:,4), 'r--');
xlabel('Time');
ylabel('Position');
title('Positions');
legend('x1', 'x2', 'y1', 'y2');
grid on;

% 2. Plot velocities (vx1, vx2, vy1, vy2)
subplot(2, 2, 2);
plot(t, X(:,5), 'b-', t, X(:,6), 'b--', t, X(:,7), 'r-', t, X(:,8), 'r--');
xlabel('Time');
ylabel('Velocity');
title('Velocities');
legend('vx1', 'vx2', 'vy1', 'vy2');
grid on;

% 3. Plot trajectories (x1 vs y1, x2 vs y2)
subplot(2, 2, 3);
plot(X(:,1), X(:,3), 'b-', X(:,2), X(:,4), 'r-');
xlabel('X Position');
ylabel('Y Position');
title('Trajectories');
legend('P1', 'P2');
grid on;
axis equal;




function dXdt = myStateSpace(t, X)
    % Extract variables from the state vector
    x_1 = X(1);
    x_2 = X(2);
    y_1 = X(3);
    y_2 = X(4);
    vx1 = X(5);
    vx2 = X(6);
    vy1 = X(7);
    vy2 = X(8);

D = 0.02;
E = 206e5;
I = (pi/64)*D^4;
rho = 1000;
Area = (pi*D^2)/4; 
mu_val = rho*Area;
L = 1;

    %Define forces (if F!=0 when calculating equations of motion)
    % if t >= 0 && t <= 0.1
    %     F = 0; % or any other function of t, e.g., sin(t)
    % else
    %     F = 0;
    % end



    %— integrand for Qben1
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



    % Define the accelerations (second derivatives)
    
    x_dotdot1 = (60*(2*Qben1 - Qben2))/pi;




    x_dotdot2 = -(20*(3*Qben1 - 4*Qben2))/pi;



    y_dotdot1 = (60*(2*Qben3 - Qben4))/pi;



    y_dotdot2 = -(20*(3*Qben3 - 4*Qben4))/pi;





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
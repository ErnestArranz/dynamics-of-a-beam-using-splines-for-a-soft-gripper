% clear
% clc
% close

%load('lastimplementation.m')

%%
% Datos iniciales
% Suponemos que q es una matriz con las coordenadas de los puntos de control
% q(:,1) -> x1, q(:,2) -> x2, q(:,3) -> y2, q(:,4) -> x3, q(:,5) -> y3, q(:,6) -> x4, q(:,7) -> y4, q(:,8) -> x5, q(:,9) -> y5
% Vector de tiempo t

% Vector de parámetro u (para la spline)
u = linspace(0, 1, 100); % Dividimos el intervalo [0,1] en 100 puntos

% Definir límites de los ejes
x_min = min([X(:,1)', X(:,2)', X(:,4)', X(:,6)', X(:,8)']); % Mínimo valor de x (incluyendo todos los puntos)
x_max = max([X(:,1)', X(:,2)', X(:,4)', X(:,6)', X(:,8)']); % Máximo valor de x (incluyendo todos los puntos)
y_min = min([X(:,3)', X(:,5)', X(:,7)', X(:,9)']);         % Mínimo valor de y (incluyendo todos los puntos)
y_max = max([X(:,3)', X(:,5)', X(:,7)', X(:,9)']);         % Máximo valor de y (incluyendo todos los puntos)

% Crear figura
figure('Position', [100, 100, 1200, 800]); % Make figure larger (width=1200, height=800)
hold on;
axis equal;
grid on;

% Inicializar trayectorias de los puntos de control
traj_P0_x = []; traj_P0_y = [];
traj_P1_x = []; traj_P1_y = [];
traj_P2_x = []; traj_P2_y = [];
traj_P3_x = []; traj_P3_y = [];
traj_P4_x = []; traj_P4_y = [];
traj_P5_x = []; traj_P5_y = []; % Añadido para el punto P5

n=100;

% Create a VideoWriter object
videoFile = 'quintic_spline_animation_force_001.avi'; % Specify the output video file name
videoWriter = VideoWriter(videoFile);
videoWriter.FrameRate = 10; % Set the frame rate (adjust as needed)
open(videoWriter);

% Animación
for k = 1:n:length(t) % skip frames to make animation faster
    % Obtener las coordenadas de los puntos de control en el tiempo t(k)
    P0 = [0, 0]; % Punto inicial (x0, y0)
    P1 = [X(k, 1), 0]; % Punto intermedio (x1, y1=0)
    P2 = [X(k, 2), X(k, 3)]; % Punto de control (x2, y2)
    P3 = [X(k, 4), X(k, 5)]; % Punto de control (x3, y3)
    P4 = [X(k, 6), X(k, 7)]; % Punto de control (x4, y4)
    P5 = [X(k, 8), X(k, 9)]; % Punto final (x5, y5) - nuevo punto
    
    % Get the current external force value

    % Guardar las trayectorias de los puntos de control
    traj_P0_x = [traj_P0_x, P0(1)]; traj_P0_y = [traj_P0_y, P0(2)];
    traj_P1_x = [traj_P1_x, P1(1)]; traj_P1_y = [traj_P1_y, P1(2)];
    traj_P2_x = [traj_P2_x, P2(1)]; traj_P2_y = [traj_P2_y, P2(2)];
    traj_P3_x = [traj_P3_x, P3(1)]; traj_P3_y = [traj_P3_y, P3(2)];
    traj_P4_x = [traj_P4_x, P4(1)]; traj_P4_y = [traj_P4_y, P4(2)];
    traj_P5_x = [traj_P5_x, P5(1)]; traj_P5_y = [traj_P5_y, P5(2)]; % Añadido para P5
    
    % Limpiar la figura anterior
    clf;
    hold on;
    axis equal;
    grid on;
    
    % Dibujar trayectorias de los puntos de control
    plot(traj_P0_x, traj_P0_y, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Trayectoria P0');
    plot(traj_P1_x, traj_P1_y, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Trayectoria P1');
    plot(traj_P2_x, traj_P2_y, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Trayectoria P2');
    plot(traj_P3_x, traj_P3_y, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Trayectoria P3');
    plot(traj_P4_x, traj_P4_y, 'c-', 'LineWidth', 1.5, 'DisplayName', 'Trayectoria P4');
    plot(traj_P5_x, traj_P5_y, 'y-', 'LineWidth', 1.5, 'DisplayName', 'Trayectoria P5'); % Añadido para P5
    
    % Dibujar puntos de control actuales
    plot(P0(1), P0(2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'Punto P0');
    plot(P1(1), P1(2), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName', 'Punto P1');
    plot(P2(1), P2(2), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', 'Punto P2');
    plot(P3(1), P3(2), 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm', 'DisplayName', 'Punto P3');
    plot(P4(1), P4(2), 'co', 'MarkerSize', 8, 'MarkerFaceColor', 'c', 'DisplayName', 'Punto P4');
    plot(P5(1), P5(2), 'yo', 'MarkerSize', 8, 'MarkerFaceColor', 'y', 'DisplayName', 'Punto P5'); % Añadido para P5

    % Dibujar líneas de control
    plot([P0(1), P1(1)], [P0(2), P1(2)], 'k:', 'LineWidth', 1);
    plot([P1(1), P2(1)], [P1(2), P2(2)], 'k:', 'LineWidth', 1);
    plot([P2(1), P3(1)], [P2(2), P3(2)], 'k:', 'LineWidth', 1);
    plot([P3(1), P4(1)], [P3(2), P4(2)], 'k:', 'LineWidth', 1);
    plot([P4(1), P5(1)], [P4(2), P5(2)], 'k:', 'LineWidth', 1); % Añadido para P5

    % Calcular y dibujar la spline completa (quíntica)
    pu_x = zeros(size(u));
    pu_y = zeros(size(u));
    for i = 1:length(u)
        % Calcular los coeficientes de Bernstein para curva quíntica
        b0 = (1 - u(i))^5;
        b1 = 5 * u(i) * (1 - u(i))^4;
        b2 = 10 * u(i)^2 * (1 - u(i))^3;
        b3 = 10 * u(i)^3 * (1 - u(i))^2;
        b4 = 5 * u(i)^4 * (1 - u(i));
        b5 = u(i)^5;
        
        % Calcular punto en la spline
        pu_x(i) = b0 * P0(1) + b1 * P1(1) + b2 * P2(1) + b3 * P3(1) + b4 * P4(1) + b5 * P5(1);
        pu_y(i) = b0 * P0(2) + b1 * P1(2) + b2 * P2(2) + b3 * P3(2) + b4 * P4(2) + b5 * P5(2);
    end
    plot(pu_x, pu_y, 'k-', 'LineWidth', 2, 'DisplayName', 'Spline');
    
    % Establecer límites de los ejes
    xlim([0, 1.2]); % Ajustar según necesidad
    ylim([-0.7, 0.7]); % Ajustar según necesidad

    % Add title with current time
    title(sprintf('t = %.2f s', t(k)));

    % Create legend with all elements
    legend({'P0 trajectory', 'P1 trajectory', 'P2 trajectory', 'P3 trajectory', 'P4 trajectory', 'P5 trajectory', ...
           ' P0', ' P1', ' P2', ' P3', ' P4', ' P5', 'Spline'}, ...
           'Location', 'northeast');
    
    % Capture the current frame and write it to the video
    frame = getframe(gcf);
    writeVideo(videoWriter, frame);

    % Pausa para crear efecto de animación
    pause(0.1);
end

% Close the VideoWriter object
close(videoWriter);

disp(['Animation saved as ', videoFile]);

%%
% Add plot comparing initial and final positions (superposed)
figure('Position', [100, 100, 1200, 800]); % Make figure larger (width=1200, height=800)
hold on;
grid on;
axis equal;

% Initial position (first time step)
P0_initial = [0, 0];
P1_initial = [X(1, 1), 0];
P2_initial = [X(1, 2), X(1, 3)];
P3_initial = [X(1, 4), X(1, 5)];
P4_initial = [X(1, 6), X(1, 7)];
P5_initial = [X(1, 8), X(1, 9)]; % Añadido para P5

% Final position (last time step)
P0_final = [0, 0];
P1_final = [X(end, 1), 0];
P2_final = [X(end, 2), X(end, 3)];
P3_final = [X(end, 4), X(end, 5)];
P4_final = [X(end, 6), X(end, 7)];
P5_final = [X(end, 8), X(end, 9)]; % Añadido para P5

% Calculate splines
u_plot = linspace(0, 1, 100);
spline_x_initial = zeros(1, length(u_plot));
spline_y_initial = zeros(1, length(u_plot));
spline_x_final = zeros(1, length(u_plot));
spline_y_final = zeros(1, length(u_plot));

for i = 1:length(u_plot)
    % Coeficientes para curva quíntica
    b0 = (1 - u_plot(i))^5;
    b1 = 5 * u_plot(i) * (1 - u_plot(i))^4;
    b2 = 10 * u_plot(i)^2 * (1 - u_plot(i))^3;
    b3 = 10 * u_plot(i)^3 * (1 - u_plot(i))^2;
    b4 = 5 * u_plot(i)^4 * (1 - u_plot(i));
    b5 = u_plot(i)^5;
    
    % Initial spline
    spline_x_initial(i) = b0 * P0_initial(1) + b1 * P1_initial(1) + b2 * P2_initial(1) + b3 * P3_initial(1) + b4 * P4_initial(1) + b5 * P5_initial(1);
    spline_y_initial(i) = b0 * P0_initial(2) + b1 * P1_initial(2) + b2 * P2_initial(2) + b3 * P3_initial(2) + b4 * P4_initial(2) + b5 * P5_initial(2);
    
    % Final spline
    spline_x_final(i) = b0 * P0_final(1) + b1 * P1_final(1) + b2 * P2_final(1) + b3 * P3_final(1) + b4 * P4_final(1) + b5 * P5_final(1);
    spline_y_final(i) = b0 * P0_final(2) + b1 * P1_final(2) + b2 * P2_final(2) + b3 * P3_final(2) + b4 * P4_final(2) + b5 * P5_final(2);
end

% Plot the initial position and spline
plot(P0_initial(1), P0_initial(2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'P0 (Initial)');
plot(P1_initial(1), P1_initial(2), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName', 'P1 (Initial)');
plot(P2_initial(1), P2_initial(2), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', 'P2 (Initial)');
plot(P3_initial(1), P3_initial(2), 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm', 'DisplayName', 'P3 (Initial)');
plot(P4_initial(1), P4_initial(2), 'co', 'MarkerSize', 8, 'MarkerFaceColor', 'c', 'DisplayName', 'P4 (Initial)');
plot(P5_initial(1), P5_initial(2), 'yo', 'MarkerSize', 8, 'MarkerFaceColor', 'y', 'DisplayName', 'P5 (Initial)');
plot(spline_x_initial, spline_y_initial, 'k--', 'LineWidth', 2, 'DisplayName', 'Initial Spline');

% Plot connecting lines for initial position
plot([P0_initial(1), P1_initial(1)], [P0_initial(2), P1_initial(2)], 'k:');
plot([P1_initial(1), P2_initial(1)], [P1_initial(2), P2_initial(2)], 'k:');
plot([P2_initial(1), P3_initial(1)], [P2_initial(2), P3_initial(2)], 'k:');
plot([P3_initial(1), P4_initial(1)], [P3_initial(2), P4_initial(2)], 'k:');
plot([P4_initial(1), P5_initial(1)], [P4_initial(2), P5_initial(2)], 'k:');

% Plot the final position and spline
plot(P0_final(1), P0_final(2), 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'P0 (Final)');
plot(P1_final(1), P1_final(2), 'gs', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName', 'P1 (Final)');
plot(P2_final(1), P2_final(2), 'bs', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', 'P2 (Final)');
plot(P3_final(1), P3_final(2), 'ms', 'MarkerSize', 8, 'MarkerFaceColor', 'm', 'DisplayName', 'P3 (Final)');
plot(P4_final(1), P4_final(2), 'cs', 'MarkerSize', 8, 'MarkerFaceColor', 'c', 'DisplayName', 'P4 (Final)');
plot(P5_final(1), P5_final(2), 'ys', 'MarkerSize', 8, 'MarkerFaceColor', 'y', 'DisplayName', 'P5 (Final)');
plot(spline_x_final, spline_y_final, 'k-', 'LineWidth', 2, 'DisplayName', 'Final Spline');

% Plot connecting lines for final position
plot([P0_final(1), P1_final(1)], [P0_final(2), P1_final(2)], 'k--');
plot([P1_final(1), P2_final(1)], [P1_final(2), P2_final(2)], 'k--');
plot([P2_final(1), P3_final(1)], [P2_final(2), P3_final(2)], 'k--');
plot([P3_final(1), P4_final(1)], [P3_final(2), P4_final(2)], 'k--');
plot([P4_final(1), P5_final(1)], [P4_final(2), P5_final(2)], 'k--');

title('Initial and Final Positions');
xlabel('X Position (m)');
ylabel('Y Position (m)');
xlim([0, 1.2]);
ylim([-0.7, 0.6]);
legend('show', 'Location', 'northeast');
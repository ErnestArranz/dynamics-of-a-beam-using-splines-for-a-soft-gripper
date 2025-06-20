%clear
clc
close
clf

%spline deformation animation code and video saver. also represents initial
%and final position. this code must be ran with the state variables in the
%workspace. that is running first dynamics.m or loading saved data

%%
% Datos iniciales
% Suponemos que q es una matriz con las coordenadas de los puntos de control
% q(:,1) -> x1, q(:,2) -> x2, q(:,3) -> y2
% Vector de tiempo t

% Vector de parámetro u (para la spline)
u = linspace(0, 1, 100); % Dividimos el intervalo [0,1] en 100 puntos

% Definir límites de los ejes
x_min = min([X(:,1)', X(:,2)']); % Mínimo valor de x (incluyendo P0, P1, P2)
x_max = max([X(:,1)', X(:,2)']); % Máximo valor de x (incluyendo P0, P1, P2)
y_min = min([X(:,3)]);         % Mínimo valor de y (incluyendo P0, P1, P2)
y_max = max([X(:,3)]);         % Máximo valor de y (incluyendo P0, P1, P2)

% Crear figura
figure;
hold on;
axis equal;
grid on;

% Inicializar trayectorias de los puntos de control
traj_P0_x = []; traj_P0_y = [];
traj_P1_x = []; traj_P1_y = [];
traj_P2_x = []; traj_P2_y = [];

n=10;

% Create a VideoWriter object
videoFile = 'spline_animation_force_001.avi'; % Specify the output video file name
videoWriter = VideoWriter(videoFile);
videoWriter.FrameRate = 10; % Set the frame rate (adjust as needed)
open(videoWriter);

% Animación
for k = 1:n:length(t) % skip frames to make animation faster%for k = 1:length(t)

%for k = 1:length(t)

%for k = 1:length(t)
    % Obtener las coordenadas de los puntos de control en el tiempo t(k)
    P0 = [0, 0]; % Punto inicial (x0, y0)
    P1 = [X(k, 1), 0]; % Punto intermedio (x1, y1=0)
    P2 = [X(k, 2), X(k, 3)]; % Punto final (x2, y2)
    
    % Get the current external force value

    % Guardar las trayectorias de los puntos de control
    traj_P0_x = [traj_P0_x, P0(1)]; traj_P0_y = [traj_P0_y, P0(2)];
    traj_P1_x = [traj_P1_x, P1(1)]; traj_P1_y = [traj_P1_y, P1(2)];
    traj_P2_x = [traj_P2_x, P2(1)]; traj_P2_y = [traj_P2_y, P2(2)];
    
    % Limpiar la figura anterior
    clf;
    hold on;
    axis equal;
    grid on;
    
    % Dibujar trayectorias de los puntos de control
    plot(traj_P0_x, traj_P0_y, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Trayectoria P0'); % Trayectoria de P0
    plot(traj_P1_x, traj_P1_y, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Trayectoria P1'); % Trayectoria de P1
    plot(traj_P2_x, traj_P2_y, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Trayectoria P2'); % Trayectoria de P2
    
    % Dibujar puntos de control actuales
    plot(P0(1), P0(2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'Punto P0'); % P0
    plot(P1(1), P1(2), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName', 'Punto P1'); % P1
    plot(P2(1), P2(2), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', 'Punto P2'); % P2

    % Calcular y dibujar la spline completa
    pu_x = zeros(size(u));
    pu_y = zeros(size(u));
    for i = 1:length(u)
        % Calcular los coeficientes de Bernstein
        b0 = (1 - u(i))^2;
        b1 = 2 * u(i) * (1 - u(i));
        b2 = u(i)^2;
        
        % Calcular punto en la spline
        pu_x(i) = b0 * P0(1) + b1 * P1(1) + b2 * P2(1);
        pu_y(i) = b0 * P0(2) + b1 * P1(2) + b2 * P2(2);
    end
    plot(pu_x, pu_y, 'k-', 'LineWidth', 2, 'DisplayName', 'Spline'); % Spline completa
    
    % Agregar leyenda
    %  legend('show', 'Location', 'northwest'); % Mostrar leyenda en la esquina superior izquierda
    
    % Establecer límites de los ejes
    xlim([0, 70]); % Agregar un margen de 1 unidad en x
    ylim([-70, 10]); % Agregar un margen de 1 unidad en y
    % ylim([-0.01, 0.01])
    %dasect[1 0.0001 1]); % Reduces compression by making y-axis larger

    % Add title with current time
    title(sprintf('t = %.2f s', t(k)));

    % Create legend with all elements
    legend({'P0 trajectory', 'P1 trajectory', 'P2 trajectory', ...
           ' P0', ' P1', ' P2', 'Spline'}, ...
           'Location', 'northwest');
    
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
figure;
hold on;
grid on;
axis equal;

% Initial position (first time step)
P0_initial = [0, 0];
P1_initial = [X(1, 1), 0];
P2_initial = [X(1, 2), X(1, 3)];

% Final position (last time step)
P0_final = [0, 0];
P1_final = [X(end, 1), 0];
P2_final = [X(end, 2), X(end, 3)];

% Calculate splines
u_plot = linspace(0, 1, 100);
spline_x_initial = zeros(1, length(u_plot));
spline_y_initial = zeros(1, length(u_plot));
spline_x_final = zeros(1, length(u_plot));
spline_y_final = zeros(1, length(u_plot));

for i = 1:length(u_plot)
    b0 = (1 - u_plot(i))^2;
    b1 = 2 * u_plot(i) * (1 - u_plot(i));
    b2 = u_plot(i)^2;
    
    % Initial spline
    spline_x_initial(i) = b0 * P0_initial(1) + b1 * P1_initial(1) + b2 * P2_initial(1);
    spline_y_initial(i) = b0 * P0_initial(2) + b1 * P1_initial(2) + b2 * P2_initial(2);
    
    % Final spline
    spline_x_final(i) = b0 * P0_final(1) + b1 * P1_final(1) + b2 * P2_final(1);
    spline_y_final(i) = b0 * P0_final(2) + b1 * P1_final(2) + b2 * P2_final(2);
end

% Plot the initial position and spline
plot(P0_initial(1), P0_initial(2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'P0 (Initial)');
plot(P1_initial(1), P1_initial(2), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName', 'P1 (Initial)');
plot(P2_initial(1), P2_initial(2), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', 'P2 (Initial)');
plot(spline_x_initial, spline_y_initial, 'k--', 'LineWidth', 2, 'DisplayName', 'Initial Spline');

% Plot connecting lines for initial position
plot([P0_initial(1), P1_initial(1)], [P0_initial(2), P1_initial(2)], 'k:');
plot([P1_initial(1), P2_initial(1)], [P1_initial(2), P2_initial(2)], 'k:');

% Plot the final position and spline
plot(P0_final(1), P0_final(2), 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'DisplayName', 'P0 (Final)');
plot(P1_final(1), P1_final(2), 'gs', 'MarkerSize', 8, 'MarkerFaceColor', 'g', 'DisplayName', 'P1 (Final)');
plot(P2_final(1), P2_final(2), 'bs', 'MarkerSize', 8, 'MarkerFaceColor', 'b', 'DisplayName', 'P2 (Final)');
plot(spline_x_final, spline_y_final, 'k-', 'LineWidth', 2, 'DisplayName', 'Final Spline');

% Plot connecting lines for final position
plot([P0_final(1), P1_final(1)], [P0_final(2), P1_final(2)], 'k--');
plot([P1_final(1), P2_final(1)], [P1_final(2), P2_final(2)], 'k--');

title('Initial and Final Positions');
xlabel('X Position (m)');
ylabel('Y Position (m)');
xlim([0 , 70]);
ylim([-70, 10]);
legend('show', 'Location', 'northwest');
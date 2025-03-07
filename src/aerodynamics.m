clc; clear; close all;
format compact;

% Kinematic Viscosity (m^2/s) at Temperature (°C)
T = [0, 15, 25]; % Temperature in Celsius
nu = [1.34E-05, 1.48E-05, 1.57E-05]; % Kinematic viscosity in m^2/s

coeffs = polyfit(T, nu, 1);

disp('Linear Fit Equation: nu(T) = ');
disp([num2str(coeffs(1), '%.5e') ' * T + ' num2str(coeffs(2), '%.5e')])

temps = [9, 18, 30] % min, avg, max, in °C
kin_visc_air = polyval(coeffs, temps) % min, avg, max, in (m^2/s)

T_fit = linspace(min([T,temps]), max([T,temps]), 100);
nu_fit = polyval(coeffs, T_fit); 

% Plot Data
figure;
hold on;
grid on;
plot(T_fit, nu_fit, 'b-', 'LineWidth', 2); % Fitted line (blue)
plot(T, nu, 'ro', 'MarkerSize', 8, 'LineWidth', 2); % Original data points (red circles)
plot(temps, kin_visc_air, 'gs', 'MarkerSize', 8, 'LineWidth', 2); % Additional points (green squares)

% Labels and Title
xlabel('Temperature (°C)');
ylabel('Kinematic Viscosity (m^2/s)');
title('Linear Fit of Kinematic Viscosity vs. Temperature');

% Legend
legend('Original Data', 'Linear Fit', 'Computed Points');

cruise_speed_names = ['eBee X Cruise', 'Switchblade Cruise', ...
    'Switchblade Dash'];
cruise_speeds = [25.2, 31.64, 51.8] % m/s

avg_chord_length_names = ['Initial Guess', 'Flying Wing', 'Giga Drone'];
avg_chord_lengths = [0.02, 0.0275, 0.04] % m

Re_values = zeros(length(cruise_speeds), length(avg_chord_lengths), length(kin_visc_air));

for i = 1:length(cruise_speeds)
    for j = 1:length(avg_chord_lengths)
        for k = 1:length(kin_visc_air)
            Re_values(i, j, k) = (cruise_speeds(i) * avg_chord_lengths(j)) / kin_visc_air(k);
        end
    end
end

% Create a meshgrid for plotting
[SpeedGrid, ViscosityGrid] = meshgrid(cruise_speeds, kin_visc_air);

for i = 1:length(avg_chord_lengths)
    % % Plotting the surface for the second chord length
    % figure;
    % surf(SpeedGrid, ViscosityGrid, squeeze(Re_values(:,:,i))); % Using the last chord length
    % xlabel('Cruise Speed (m/s)');
    % ylabel('Kinematic Viscosity (m^2/s)');
    % zlabel('Reynolds Number');
    % title(['Reynolds Number for Chord Length: ' num2str(avg_chord_lengths(i)) ' meters']);
    % 
    % % Display the plot
    % shading interp;  % Smooth the colors of the surface
    % colorbar;  % Add color bar to indicate the Reynolds number magnitude
    % grid on;

    % Create a meshgrid for plotting
    [SpeedGrid, TempGrid] = meshgrid(cruise_speeds, temps);
    
    % Plotting the surface for the second chord length
    figure;
    surf(SpeedGrid, TempGrid, squeeze(Re_values(:,:,i))); % Using the second chord length
    xlabel('Cruise Speed (m/s)');
    ylabel('Temperature (°C)');  % Label temperature on the Y-axis
    zlabel('Reynolds Number');
    title(['Reynolds Number for Chord Length: ' num2str(avg_chord_lengths(i)) ' meters']);
    
    % Display the plot
    shading interp;  % Smooth the colors of the surface
    colorbar;  % Add color bar to indicate the Reynolds number magnitude
    grid on;
end
% Reynolds Number
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

% Calculate reynolds numbers for each cruise speed,
% chord length, and kin visc of air.
for i = 1:length(cruise_speeds)
    for j = 1:length(avg_chord_lengths)
        for k = 1:length(kin_visc_air)
            Re_values(i, j, k) = (cruise_speeds(i) * avg_chord_lengths(j)) / kin_visc_air(k);
            % fprintf('Cruise Speed: %.2f\nChord Length: %.2f\nKin Visc Air: %.2f\nReynolds Number: %.2f\n\n', cruise_speeds(i), avg_chord_lengths(j), kin_visc_air(k), Re_values(i,j,k))
        end
    end
end

% Plot Reynolds 3D maps keeping chord length constant
% for j = 1:length(avg_chord_lengths)
%     figure;
% 
%     % Generate meshgrid for cruise speed and temperature
%     [SpeedGrid, TempGrid] = meshgrid(cruise_speeds, temps);
% 
%     % Extract correct slice of Reynolds number matrix
%     Re_surface = squeeze(Re_values(:, j, :))';  % Transpose to align correctly
% 
%     % Surface plot
%     surf(SpeedGrid, TempGrid, Re_surface);
% 
%     xlabel('Cruise Speed (m/s)');
%     ylabel('Temperature (°C)');
%     zlabel('Reynolds Number');
%     title(['Reynolds Number for Chord Length: ' num2str(avg_chord_lengths(j)) ' m']);
% 
%     % Visualization improvements
%     shading interp;
%     colorbar;
%     grid on;
% end
%% Forces
clc; clear; close all;
format compact;

% LIFT CALCULATIONS
% Air Density (kg/m^3) at Temperature (°C)
T = [0, 15, 25, 30]; % Temperature in Celsius
rho_data = [1.292, 1.225, 1.184, 1.164]; % Density in [kg/m^3]

% Perform Linear Fit
coeffs = polyfit(T, rho_data, 1);

% Display Linear Fit Equation
% disp('Linear Fit Equation: rho(T) = ');
% disp([num2str(coeffs(1), '%.5e') ' * T + ' num2str(coeffs(2), '%.5e')]);

% Compute Density for Given Temperatures
temps = [9, 18];
rho_data_2 = polyval(coeffs, temps);

% Generate Data for the Fit Line
T_fit = linspace(min(T), max(T), 100);
rho_fit = polyval(coeffs, T_fit);

% % Plot Data
% figure;
% hold on;
% grid on;
% plot(T_fit, rho_fit, 'b-', 'LineWidth', 2); % Fitted line (blue)
% plot(T, rho_data, 'ro', 'MarkerSize', 8, 'LineWidth', 2); % Original data points (red circles)
% plot(temps, rho_data_2, 'gs', 'MarkerSize', 8, 'LineWidth', 2); % Additional points (green squares)
% 
% % Labels and Title
% xlabel('Temperature (°C)');
% ylabel('Air Density (kg/m^3)');
% title('Linear Fit of Air Density vs. Temperature');
% 
% % Legend (Reordered for clarity)
% legend('Linear Fit', 'Original Data', 'Computed Points', 'Location', 'SouthWest');

% Airfoil: http://airfoiltools.com/airfoil/details?airfoil=df102-il
Cl = 0.25; % [dimensionless], measured at alpha = 0.005
V = 25.2; % [m/s], chosen from Reynolds 3D maps

% Density of Air at 18°C
T = 18;
rho = polyval(coeffs, T); % [kg/m3]
fprintf('Density of Air at 18°C:\nT = %.2f [°C]\nρ = %.2f [kg/m3]\n\n', T, rho)

% Average Wing Chord Length
c = 0.4; % [m]
% Wing span
b = 1.3; % [m]
% Wing area
A = c*b; % [m^2]

% Lift
FL = 1/2*rho*V^2*A*Cl; % N

fprintf(['Average Wing Chord Length\n' ...
    '  c = %.2f [m]\n' ...
    'Wing Span:\n' ...
    '  b = %.2f [m]\n' ...
    'Wing Area:\n' ...
    '  A = %.2f [m^2]\n' ...
    'Coefficient of Lift:\n' ...
    '  Cl = %.2f\n' ...
    'Velocity:\n' ...
    '  V = %.2f\n' ...
    'Lift Force:\n' ...
    '  L = 1/2*rho*V^2*A*Cl = %.2f [N]\n\n'], c, b, A, Cl, V, FL)

% Weight of drone
m_max = 1.5; % [kg]

m_payload = 0.4; % [kg]
m_motor = 0.175; % [kg]
m_servo = 0.018; % [kg]
m_battery = 0.152; % [kg]

m_comp = m_payload + m_motor + 3*m_servo + m_battery;
m_body = m_max - m_comp;
m = m_comp + m_body;

g = 9.82; % m/s^2

FW = m*g;

fprintf(['Weight Calculations\n' ...
    'Max Weight:\n' ...
    '  m_max = %.2f [kg]\n' ...
    'Payload Mass:\n' ...
    '  m_payload = %.2f [kg]\n' ...
    'Motor Mass:\n' ...
    '  m_motor = %.2f [kg]\n' ...
    'Servo Mass:\n' ...
    '  m_servo = %.2f [kg]\n' ...
    'Battery Mass:\n' ...
    '  m_battery = %.2f [kg]\n\n' ...
    'Mass of Components:\n' ...
    '  m_comp = m_payload + m_motor + 3*m_servo + m_battery = %.2f [kg]\n' ...
    'Mass of Body:\n' ...
    '  m_body = m_max - m_comp = %.2f [kg]\n' ...
    'Total Weight:\n' ...
    '  m = m_max = m_comp + m_body = %.2f [kg]\n\n' ...
    'Gravity:\n' ...
    '  g = %.2f [m/s^2]\n' ...
    'Weight:\n' ...
    '  W = m*g = %.2f [N]\n\n'], ...
    m_max, m_payload, m_motor, m_servo, m_battery, m_comp, m_body, m, g, FW)

% Lift-to-weight ratio
LWR = FL/FW;

% Is the lift enough to counteract weight of drone?
if (FL>FW)
    fprintf('Yes, we have enough lift!')
else
    fprintf('No, we do not have enough lift!')
end
fprintf('\nL/W Ratio: %.2f\n\n', LWR)

% DRAG CALCULATIONS

Cd = 0.01; % Extracted from Cl/Cd graph.
% Drag
FD = 1/2*rho*V^2*A*Cd % N

fprintf(['Drag Calculations:\n' ...
    'Coefficient of Drag:\n' ...
    '  Cd = %.2f\n' ...
    'Drag:\n' ...
    '  D = 1/2*rho*V^2*A*Cd = %.2f [N]\n\n'], ...
    Cd, FD)

% THRUST CALCULATIONS
m = 2.3; % Max "thrust" of motor
% Static thrust
T_static = m * g; % N, when static. Theoretical number, will be smaller when in flight.

% Thrust-to-drag ratio
TDR = T/FD;

fprintf(['Thrust Calculations:\n' ...
    'Max Thrust of Motor:\n' ...
    '  m = %.2f\n' ...
    'Static Thrust:\n' ...
    '  T_static = m*g = %.2f\n\n'], ...
    m, T_static)

% Is the thrust enough to counteract drag of drone?
if (T>FD)
    disp('Yes, we have enough thrust!')
else
    disp('No, we do not have enough thrust!')
end
fprintf('T/D Ratio: %.2f\n', TDR)

% Export Variables
save("aerodynamics.mat")
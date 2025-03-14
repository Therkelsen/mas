% Reynolds Number
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

T = 18;
rho = polyval(coeffs, T); % Density of Air at 18°C [kg/m3]
fprintf('Density of Air at %i°C:\n  ρ = %.2f [kg/m3]\n', T, rho)

b = 1; % Wing span [m]
c = 0.25; % Chord Length [m]
A = c*b; % Wing Area [m^2]
fprintf('\nWing Area = Wing Span x Chord Length\n  A = b x c\n    = %i [m] x %.2f [m]\n    = %.2f [m^2]\n', b, c, A)


% Weight of drone
g = 9.82; % Gravity [N]

m_max = 1.5; % [kg]

m_payload = 0.4; % [kg]
m_motor = 0.175; % [kg]
m_servo = 0.018; % [kg]
m_battery = 0.152; % [kg]

m_comp = m_payload + m_motor + 3*m_servo + m_battery;
m_body = m_max - m_comp;
m = m_comp + m_body;

FW = m*g;

fprintf(['\nWEIGHT STUFF\nGravity:\n  g = %.2f [N]\n' ...
    'Weight:\n  m = %.2f [kg]\n' ...
    '\nWeight Force:\n  FW = m x g\n' ...
    '     = %.2f x %.2f\n' ...
    '     = %.2f [N]\n'], g, m, m, g, FW)

% STALL SPEED
Clstall = 1.35; % Coefficient of Lift at Stall
Cdstall = 0.05; % Coefficient of Drag at Stall
alphastall = 14; % Degrees

% Stall Speed
Vstall = sqrt((2*m*g)/(rho*A*Clstall));

fprintf(['\nLIFT STUFF\nLift Coefficient at Stall:\n  Clstall = %.2f\n' ...
    'Drag Coefficient at Stall:\n  Cldrag = %.2f\n'], Clstall, Cdstall)

% % Lift
% Airfoil: http://airfoiltools.com/airfoil/details?airfoil=df102-il
Cl = 0.4; % [dimensionless], measured at alpha = 0%
% FL = 1/2*rho*V^2*A*Cl; % N
FL = FW;
Vcruise = sqrt((2*FL)/(rho*A*Cl));

fprintf(['\nLift Coefficient:\n  Cl = %.2f\n' ...
    'Lift Force = Weight Force\n' ...
    '  FL = FW = %.2f [N]\n' ...
    'Cruise Speed:\n  Vcruise = %.2f [m/s]\n'], Cl, FL, Vcruise)

% DRAG CALCULATIONS
Cd = 0.01; % Extracted from Cl/Cd graph.
% Drag
FD = 1/2*rho*Vcruise^2*A*Cd; % N

fprintf(['\nDRAG STUFF\nDrag Coefficient:\n  Cd = %.2f\n' ...
    'Drag Force:\n  FD = %.2f [N]\n'], Cd, FD)

% THRUST CALCULATIONS
m = 2.3; % Max "thrust" of motor
% Static thrust
Tstatic = m * g % N, when static. Theoretical number, will be smaller when in flight.
eta = 0.8; % Efficiency
Tactual = Tstatic * eta

% Thrust-to-drag ratio
TDR = Tactual/FD;

fprintf(['\nTHRUST STUFF\nMax Thrust Spec:\n  m = %.2f' ...
    '\nStatic Thrust:\n  Tstatic = m x g\n          = %.2f x %.2f\n          = %.2f [N]' ...
    '\nEfficiency:\n  η = %.2f = %i%%' ...
    '\nActual Thrust:\n  Tactual = Tstatic x η\n          = %.2f [N]\n'], m, m, g, Tstatic, eta, eta*100, Tactual)

% Is the thrust enough to counteract drag of drone?
if (Tactual>FD)
    fprintf('\nYes, we have enough thrust!')
else
    fprintf('\nNo, we do not have enough thrust!')
end
fprintf('\nT/D Ratio: %.2f\n', TDR)

% % REYNOLDS STUFF
% Kinematic Viscosity (m^2/s) at Temperature (°C)
T = [0, 15, 25]; % Temperature in Celsius
nu = [1.34E-05, 1.48E-05, 1.57E-05]; % Kinematic viscosity in m^2/s

coeffs = polyfit(T, nu, 1);

% disp('Linear Fit Equation: nu(T) = ');
% disp([num2str(coeffs(1), '%.5e') ' * T + ' num2str(coeffs(2), '%.5e')])

temps = [9, 18, 30]; % min, avg, max, in °C
kin_visc_air = polyval(coeffs, temps); % min, avg, max, in (m^2/s)

T_fit = linspace(min([T,temps]), max([T,temps]), 100);
nu_fit = polyval(coeffs, T_fit);

T = temps(2);
nu = nu_fit(2);

fprintf(['\nKinematic Viscosity at T = %i °C:\n  ν = %i [m^2/s]'], T, nu)

Re = (Vcruise * c) / nu;

fprintf('\nReynolds Number:\n  Re = %i', Re)

Pmotor = (Tstatic * Vcruise) / eta; % Motor Power [W]

PWR = Pmotor / m; % Power-to-weight ratio

Vmax = sqrt((2*Tactual)/(rho*Cd*A));

% Export Variables
save("aerodynamics.mat")
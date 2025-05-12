%% Aerodynamics for Steady Level Flight (All velocities are constant)
clc; clear; close all;
format compact;

%% Reynolds Stuff
% Kinematic Viscosity (ν) in m^2/s at corresponding Temperature (°C)
temp_data = [0, 15, 25]; 
nu_data = [1.34E-05, 1.48E-05, 1.57E-05];

% Fit linear model to viscosity data
coeffs = polyfit(temp_data, nu_data, 1); 

% Estimate kinematic viscosity at given temperature conditions
T_range = [9, 18, 30]; % [min, avg, max] temperature in °C
kin_visc_air = polyval(coeffs, T_range); 

% Cruise condition parameters
T_cruise = T_range(2);         % Avg temperature
nu_cruise = kin_visc_air(2);   % ν at average temp

% Generate smooth line for fitted model
T_fit = linspace(min([temp_data, T_range])-5, max([temp_data, T_range])+5, 100);
nu_fit = polyval(coeffs, T_fit);

% Plot
figure;
plot(T_fit, nu_fit, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Linear fit');
hold on;
plot(temp_data, nu_data, 'o', 'MarkerFaceColor', 'b', 'DisplayName', 'Measured data');
plot(T_range, kin_visc_air, 'ks', 'MarkerFaceColor', 'g', 'DisplayName', 'Estimated values');

xlabel('Temperature (°C)');
ylabel('Kinematic Viscosity (m^2/s)');
title('Kinematic Viscosity of Air vs Temperature');
legend('Location', 'northwest');
grid on;

b = 0.65;       % Wing span [m]
c = 0.25;       % Chord Length [m]
t = 0.025;      % Max thickness of airfoil [m]
% h_spar = 0.025; % Spar height [m]
h_spar = (0.032 + 0.026)/2; % Spar height [m]
Awing = c * b;  % Wing Area [m^2]
Vwing = Awing * t; % Wing Volume [m^3]
vcruise = 18;   % Cruise Speed [m/s]

% Calculate Reynolds number
Re = (vcruise * c) / nu_cruise;

% Display
fprintf(['\nWing Area = Wing Span x Chord Length\n' ...
    '  A = b x c\n    = %.2f [m] x %.2f [m]\n' ...
    '    = %.4f [m^2]\n' ...
    'Wing Volume = Wing Area x Wing Thickness\n' ...
    '  V = A x t\n    = %.4f [m^2] x %.3f [m]\n' ...
    '    = %.6f [m^3]\n'], ...
    b, c, Awing, Awing, t, Vwing);

fprintf(['\nKinematic Viscosity at T = %.1f °C:\n  ν = %.2e [m^2/s]' ...
    '\nCruise Speed:\n  vcruise = %.2f [m/s]\n' ...
    '\nReynolds Number:\n  Re = %.3e\n'], ...
    T_cruise, nu_cruise, vcruise, Re);

%% Weight Stuff
g = 9.82; % Gravity [N]

m_max = 1.5; % [kg]

m_payload = 0.4; % [kg]
m_motor = 0.175; % [kg]
m_servo = 0.024; % [kg]
% m_battery = 0.074; % [kg]
m_battery = 0.180*2;

m_comp = m_payload + m_motor + 3*m_servo + m_battery;
m_body = m_max - m_comp;
m = m_comp + m_body;

FW = m*g;

fprintf(['\nWEIGHT\nGravity:\n  g = %.2f [N]\n' ...
    'Weight:\n  m = %.2f [kg]\n' ...
    '\nWeight Force:\n  FW = m x g\n' ...
    '     = %.2f x %.2f\n' ...
    '     = %.2f [N]\n'], g, m, m, g, FW)

%% Lift Stuff
T = 18; % Temperature [°C]
h = 120; % Altitude [m]
[rho, P] = air_density(T, h); % Density of Air [kg/m^3]
fprintf(['\nLIFT\nDensity of Air at %i°C:\n  ρ = %.4f [kg/m3]\n' ...
    'Pressure of Air at %i°C:\n  P = %.4f [hPa]\n'], T, rho, T, P/100)

% Airfoil: http://airfoiltools.com/airfoil/details?airfoil=sd7037-il
% Cl = 0.4; % [dimensionless], measured at alpha = 0%
% FL = 1/2*rho*vcruise^2*A*Cl; % N
FL = FW;
% vcruise = sqrt((2*FL)/(rho*A*Cl));
Cl = (2*FL)/(rho*vcruise^2*Awing);

fprintf(['Lift Force:\n' ...
    '  FL = %.2f [N]\n' ...
    'Lift Coefficient:\n  Cl = %.2f\n' ...
    ], FL, Cl)

%% Airfoil Stuff
airfoil = readtable('xf-sd7037-il-500000.csv', 'HeaderLines', 10);

% Target lift coefficient
Cd = 0; % Drag Coefficient
alpha = 0; % Angle of Attack [°]
tol = 0.01;

% Find matching index
indices = find(abs(airfoil.Cl - Cl) <= tol);

fprintf('\nAIRFOIL')
if ~isempty(indices)
    Cd = airfoil.Cd(indices(1));
    alpha = airfoil.Alpha(indices(1));
    fprintf(['\nFor Cl = %.4f,\n' ...
        'Cd = %.4f\n' ...
        'alpha = %.1f [°]\n'], Cl, Cd, alpha)
else
    fprintf('\nNo matching Cl found within tolerance.\n');
end

% Extract relevant data
Cl_data = airfoil.Cl;
Cd_data = airfoil.Cd;
alpha_data = airfoil.Alpha;

% Plot Cl vs Cd
figure;
plot(Cd_data, Cl_data, 'b.-', 'HandleVisibility', 'off');
hold on;
if ~isempty(indices)
    xline(Cd, 'k--', 'HandleVisibility', 'off');
    yline(Cl, 'k--', 'HandleVisibility', 'off');
    plot(Cd, Cl, 'ks', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', sprintf('(C_L, C_D) = (%.2f, %.2f)', Cl, Cd));
end
xlabel('C_D');
ylabel('C_L');
title('Lift Coefficient vs Drag Coefficient');
legend('Location', 'best');
grid on;

% Plot Cl vs Alpha
figure;
plot(alpha_data, Cl_data, 'r.-', 'HandleVisibility', 'off');
hold on;
if ~isempty(indices)
    xline(alpha, 'k--', 'HandleVisibility', 'off');
    yline(Cl, 'k--', 'HandleVisibility', 'off');
    plot(alpha, Cl, 'ks', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', sprintf('(C_L, \\alpha) = (%.2f, %.2f°)', Cl, alpha));
end
xlabel('\alpha [°]');
ylabel('C_L');
title('Lift Coefficient vs Angle of Attack');
legend('Location', 'best');
grid on;

%% Drag Stuff
FD = 1/2*rho*vcruise^2*Awing*Cd; % Drag Force [N]

fprintf(['\nDRAG\nDrag Force:\n  FD = %.2f [N]\n'], FD)

%% Thrust Calculation
motor = readtable('motor_data.csv', 'TextType', 'string');

idx = strcmp(motor.Parameter, "Static_Thrust");
T = motor.Value(idx);
FT = T*g/9.82/1000;

fprintf(['\nTHRUST\nThrust Force:\n  FT = %.2f [N]\n' ...
    'Thrust Weight:\n  T = %.2f [kg]\n'], FT, T)

% Thrust to Drag ratio
TDR = FT/FD;

% Is our thrust higher than our drag?
if (FT>FD)
    fprintf('\nYes, our thrust is higher than our drag!')
else
    fprintf('\nNo, our thrust is not higher than our drag!')
end
fprintf('\nThrust/Drag Ratio: %.2f\n\n', TDR)

%% Stall Stuff
Clstall = 1.3934; % Coefficient of Lift at Stall
Cdstall = 0.02657; % Coefficient of Drag at Stall
alphastall = 11.5; % Degrees

% m_actual = 1.430;
m_actual = 1.615;

m = m_actual;

% Stall Speed
vstall = sqrt((2*m*g)/(rho*Awing*Clstall));

fprintf(['\nSTALL\nLift Coefficient at Stall:\n  Clstall = %.2f\n' ...
    'Drag Coefficient at Stall:\n  Cldrag = %.2f\n' ...
    'Stall Speed: \n  vstall = %.2f\n'], Clstall, Cdstall, vstall)

% Cruise to Stall ratio
CSR = vcruise/vstall;

% Is our cruise speed higher than our stall speed?
if (vcruise>vstall)
    fprintf('\nYes, our cruise speed is\nhigh enough to not stall!')
else
    fprintf('\nNo, our cruise speed is\nnot high enough to not stall!')
end
fprintf('\nCruise/Stall Ratio: %.2f\n\n', CSR)

%% More Thrust
% Aircraft performance estimation script

% Inputs
cell_voltage_nominal = 3.7;
cell_voltage_charged = 4.2;
num_cells_series = 3;

battery_capacity_mAh = 2200;      % mAh
battery_C = 30;                   % C-rating

motor_Kv = 870;                   % RPM/V
motor_max_current = 41;           % Amps
motor_thrust_ref = 2300;          % grams at 8700 RPM
motor_RPM_ref = 8700;             % RPM at ~10V
motor_efficiency = 0.8;           % 0-1

h_body = 0.08; % [m] outer diameter
l_body = 0.45; % [m] length of body

esc_max_current = 40;

prop_diameter_in = 12;            % inches
prop_pitch_in = 6;                % inches

Cd_plane = 0.08;

static_to_actual_thrust = 0.65;

% Air density
rho = 1.225;  % kg/m³ at sea level

% Calculated battery voltages
battery_voltage_nominal = cell_voltage_nominal * num_cells_series;
battery_voltage_charged = cell_voltage_charged * num_cells_series;

% Motor RPM and thrust estimation
motor_RPM = motor_Kv * battery_voltage_charged;

% Adjust the thrust scaling (use a conservative exponent)
alpha = 2;  % exponent to adjust thrust scaling

scaled_thrust_g = motor_thrust_ref * (motor_RPM / motor_RPM_ref)^alpha;
scaled_thrust_N = scaled_thrust_g / 1000 * g * static_to_actual_thrust;  % convert thrust to Newtons

% Current estimation (linear scaling with thrust)
motor_current = motor_max_current * (scaled_thrust_g / motor_thrust_ref);
motor_current = min(motor_current, esc_max_current);   % ESC current limit

% Battery current capability and estimated flight time
battery_max_current = (battery_capacity_mAh / 1000) * battery_C;
flight_time_minutes = (battery_capacity_mAh / 1000) / motor_current * 60;

% Propeller pitch speed estimation (ideal airspeed)
pitch_speed_mps = (motor_RPM * prop_pitch_in * 0.0254) / 60;

kmh_to_mps = 0.277778;

propcalc_cruise_speed = 94*870/985*kmh_to_mps; % From propcalc
propcalc_rpm = 9967;
propcalc_pitch_speed = 81*kmh_to_mps;
pitch = propcalc_pitch_speed*60/propcalc_rpm;
K = propcalc_rpm*pitch/propcalc_cruise_speed;
cruise_speed = (motor_RPM * pitch)/K;

% --- Print Results ---
fprintf('\n--- Battery ---\n');
fprintf('Nominal Voltage: %.2f V\n', battery_voltage_nominal);
fprintf('Charged Voltage: %.2f V\n', battery_voltage_charged);
fprintf('Max Battery Current (%.0fC): %.1f A\n', battery_C, battery_max_current);

fprintf('\n--- Motor ---\n');
fprintf('Estimated RPM: %.0f RPM\n', motor_RPM);
fprintf('Estimated Current Draw: %.1f A\n', motor_current);
fprintf('Estimated Thrust: %.1f g (%.2f N)\n', scaled_thrust_g, scaled_thrust_N);

fprintf('\n--- Propeller ---\n');
fprintf('Estimated Pitch Speed: %.1f m/s\n', pitch_speed_mps);

fprintf('\n--- Flight Performance ---\n');
fprintf('Estimated Flight Time at Full Power: %.1f minutes\n', flight_time_minutes);
fprintf('Estimated Cruise Speed: %.2f m/s\n', cruise_speed);

%% Aileron Sizing
% === Aileron Deflection and Position ===
delta_a = deg2rad(20)      % Aileron deflection angle [rad]
aileron_span_ratio = 0.25;  % Ratio of wing span covered by each aileron

% === Wing Geometry ===
y_inner = b/2 * (1 - aileron_span_ratio)
y_outer = b/2
y_avg = (y_inner + y_outer) / 2   % Average distance from centerline

% === Moment of Inertia ===
Ix = (1/12) * m * (b^2); % Estimate of moment of inertia about x-axis [kg·m^2]

% === Desired Roll Rate ===
desired_roll_rate = deg2rad(200)  % Desired roll rate [rad/s]

% === Roll Moment Required ===
% Roll moment needed to reach the desired roll rate
L_required = desired_roll_rate * Ix;

% === Solve for Required Aileron Area ===
% Rearranged formula for Cl_delta:
% Cl = L / (0.5 * rho * vcruise^2 * Awing * b)
Cl_required = L_required / (0.5 * rho * vcruise^2 * Awing * b);

% Cl_delta approximation:
% Cl_delta ≈ (pi/4) * (Sa * y_avg) / (Awing * b/2)
% Solve for Sa (aileron area)
Sa = (Cl_required / delta_a) * (4 / pi) * (Awing * b/2) / y_avg;

% === Convert Aileron Area to Dimensions ===
% Assume aileron spans aileron_span_ratio of half the wing span
aileron_span = b/2 * aileron_span_ratio;
aileron_chord = Sa / aileron_span;

% === Output Results ===
fprintf('AILERON SIZING\n');
fprintf('Required roll moment: %.2f Nm\n', L_required);
fprintf('Required Cl: %.4f\n', Cl_required);
fprintf('Estimated aileron area: %.4f m^2\n', Sa);
fprintf('Suggested aileron span: %.3f m\n', aileron_span);
fprintf('Suggested aileron chord: %.3f m\n', aileron_chord);
fprintf('Deflection used: %.1f deg\n', rad2deg(delta_a));


%% Export Variables
save("aerodynamics.mat")

%% Functions
function [rho, P] = air_density(T_celsius, h)
    % Constants
    P0 = 101325; % Sea level standard atmospheric pressure (Pa)
    T0 = 288.15; % Sea level standard temperature (K)
    g = 9.82; % Gravitational acceleration (m/s²)
    M = 0.0289644; % Molar mass of Earth's air (kg/mol)
    R_univ = 8.3144598; % Universal gas constant (J/(mol·K))
    R_air = 287.05; % Specific gas constant for dry air (J/(kg·K))
    L = 0.0065; % Temperature lapse rate (K/m)
    
    % Convert temperature to Kelvin
    T = T_celsius + 273.15;
    
    % Calculate pressure at altitude using barometric formula
    P = P0 * (1 - (L * h) / T0)^( (g * M) / (R_univ * L) );

    % Calculate air density
    rho = P / (R_air * T);
end
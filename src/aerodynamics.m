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

b = 0.65;       % Wing span [m]
c = 0.25;       % Chord Length [m]
t = 0.025;      % Max thickness of airfoil [m]
Awing = c * b;  % Wing Area [m^2]
Vwing = Awing * t; % Wing Volume [m^3]
vcruise = 25;   % Cruise Speed [m/s]

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
m_servo = 0.018; % [kg]
m_battery = 0.074; % [kg]

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
Cd = 0; % Drag Coefficient
alpha = 0; % Angle of Attack [°]
tol = 0.005;
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
Clstall = 1.35; % Coefficient of Lift at Stall
Cdstall = 0.05; % Coefficient of Drag at Stall
alphastall = 14; % Degrees

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

%% Aileron Sizing
Ix = (1/12) * m * (b^2); % Estimate of moment of inertia about x-axis [kg·m^2]

desired_roll_rate = deg2rad(200);  % Desired roll rate [rad/s]

% Aileron deflection angle (max)
delta_a = deg2rad(20);      % Aileron deflection angle [rad]

% Assumed aileron position
aileron_span_ratio = 0.25;  % Ratio of wing span covered by each aileron
y_inner = b/2 * (1 - aileron_span_ratio);
y_outer = b/2;
y_avg = (y_inner + y_outer) / 2;   % Average distance from centerline

% === Roll Moment Required ===
% Roll moment needed to reach the desired roll rate
L_required = desired_roll_rate * Ix;

% === Solve for required aileron area ===
% Rearranged formula for Cl_delta:
% Cl = L / (0.5 * rho * vcruise^2 * Awing * b)
Cl_required = L_required / (0.5 * rho * vcruise^2 * Awing * b);

% Cl_delta approximation:
% Cl_delta ≈ (pi/4) * (Sa * y_avg) / (Awing * b/2)
% Solve for Sa (aileron area)
Sa = (Cl_required / delta_a) * (4 / pi) * (Awing * b/2) / y_avg;

% === Convert area to dimensions ===
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
    g = 9.80665; % Gravitational acceleration (m/s²)
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
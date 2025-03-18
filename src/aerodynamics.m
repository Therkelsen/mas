%%
clc; clear; close all;
format compact;

%% Weight Stuff
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

fprintf(['WEIGHT\nGravity:\n  g = %.2f [N]\n' ...
    'Weight:\n  m = %.2f [kg]\n' ...
    '\nWeight Force:\n  FW = m x g\n' ...
    '     = %.2f x %.2f\n' ...
    '     = %.2f [N]\n'], g, m, m, g, FW)

%% Lift Stuff
T = 18; % Temperature [°C]
h = 120; % Altitude [m]
rho = air_density(T, h); % Density of Air [kg/m^3]
fprintf('\nLIFT\nDensity of Air at %i°C:\n  ρ = %.4f [kg/m3]\n', T, rho)

b = 1; % Wing span [m]
c = 0.25; % Chord Length [m]
A = c*b; % Wing Area [m^2]
fprintf('\nWing Area = Wing Span x Chord Length\n  A = b x c\n    = %i [m] x %.2f [m]\n    = %.2f [m^2]\n', b, c, A)

% Airfoil: http://airfoiltools.com/airfoil/details?airfoil=sd7037-il
% Cl = 0.4; % [dimensionless], measured at alpha = 0%
Vcruise = 25; % Cruise Speed [m/s]
% FL = 1/2*rho*Vcruise^2*A*Cl; % N
FL = FW;
% Vcruise = sqrt((2*FL)/(rho*A*Cl));
Cl = (2*FL)/(rho*Vcruise^2*A);

fprintf(['\nCruise Speed:\n  Vcruise = %.2f [m/s]\n' ...
    'Lift Force:\n' ...
    '  FL = %.2f [N]\n' ...
    'Lift Coefficient:\n  Cl = %.2f\n' ...
    ], Vcruise, FL, Cl)

%% Reynolds Stuff
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

fprintf('\nReynolds Number:\n  Re = %i\n', Re)

%% Airfoil Stuff
airfoil = readtable('xf-sd7037-il-500000.csv', 'HeaderLines', 10);
Cd = 0; % Drag Coefficient
alpha = 0; % Angle of Attack [Degrees]
tol = 0.005;
indices = find(abs(airfoil.Cl - Cl) <= tol);
fprintf('\nAIRFOIL')
if ~isempty(indices)
    Cd = airfoil.Cd(indices(1));
    alpha = airfoil.Alpha(indices(1));
    fprintf(['\nFor Cl = %.4f,\n' ...
        'Cd = %.4f\n' ...
        'alpha = %.1f\n'], Cl, Cd, alpha)
else
    fprintf('\nNo matching Cl found within tolerance.\n');
end

%% Drag Stuff
FD = 1/2*rho*Vcruise^2*A*Cd; % Drag Force [N]

fprintf(['\nDRAG\nDrag Force:\n  FD = %.2f [N]\n'], FD)

%% Thrust Stuff
m = 2.3; % Max "thrust" of motor
% Static thrust
Tstatic = m * g; % N, when static. Theoretical number, will be smaller when in flight.
eta = 0.8; % Efficiency
Tactual = Tstatic * eta;

% Thrust-to-drag ratio
TDR = Tactual/FD;

fprintf(['\nTHRUST\nMax Thrust Spec:\n  m = %.2f' ...
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

%% Stall Stuff
Clstall = 1.35; % Coefficient of Lift at Stall
Cdstall = 0.05; % Coefficient of Drag at Stall
alphastall = 14; % Degrees

% Stall Speed
Vstall = sqrt((2*m*g)/(rho*A*Clstall));

fprintf(['\nSTALL\nLift Coefficient at Stall:\n  Clstall = %.2f\n' ...
    'Drag Coefficient at Stall:\n  Cldrag = %.2f\n' ...
    'Stall Speed: \n  Vstall = %.2f\n'], Clstall, Cdstall, Vstall)

%% Max Speed Stuff
Pmotor = (Tstatic * Vcruise) / eta; % Motor Power [W]

PWR = Pmotor / m; % Power-to-weight ratio

Vmax = sqrt((2*Tactual)/(rho*Cd*A));

%% Export Variables
save("aerodynamics.mat")

%% Functions
function rho = air_density(T_celsius, h)
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
    
    % Calculate pressure at altitude using barometric formula (valid for h < 11 km)
    P = P0 * (1 - (L * h) / T0)^( (g * M) / (R_univ * L) );

    % Calculate air density using the ideal gas law
    rho = P / (R_air * T);
end
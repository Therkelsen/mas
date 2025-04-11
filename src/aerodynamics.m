%% Aerodynamics for Steady Level Flight (All velocities are constant)
clc; clear; close all;
format compact;

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

b = 0.65; % Wing span [m]
c = 0.25; % Chord Length [m]
Awing = c*b; % Wing Area [m^2]
t = 0.025; % Airfoil Maximum Thickness [m]
vcruise = 25; % Cruise Speed [m/s]
fprintf(['\nWing Area = Wing Span x Chord Length\n' ...
    '  A = b x c\n    = %i [m] x %.2f [m]\n' ...
    '    = %.2f [m^2]\n'], ...
    b, c, Awing)

Re = (vcruise * c) / nu;

fprintf(['\nKinematic Viscosity at T = %i °C:\n  ν = %i [m^2/s]' ...
    '\nCruise Speed:\n  vcruise = %.2f [m/s]\n' ...
    '\nReynolds Number:\n  Re = %i\n'], T, nu, vcruise, Re)

%% Weight Stuff
g = 9.82; % Gravity [N]

m_max = 1.5; % [kg]

m_payload = 0.4; % [kg]
m_motor = 0.175; % [kg]
m_servo = 0.018; % [kg]
m_battery = 0.120; % [kg]

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
rho = air_density(T, h); % Density of Air [kg/m^3]
fprintf('\nLIFT\nDensity of Air at %i°C:\n  ρ = %.4f [kg/m3]\n', T, rho)

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
Dprop = 12 * 0.0254; % Propeller diameter [m] (12 in converted to m)
rprop = Dprop / 2; % Propeller radius [m]
Aprop = pi * rprop^2; % Propeller disk area [m^2]

FTmax = 2300 * g / 1000; % Convert max thrust from grams to Newtons

% Compute exit velocity using momentum theory
ve = sqrt(2 * FTmax / (rho * Aprop)); % Exit velocity [m/s]

% Compute thrust using the momentum equation
FT = 1/2 * rho * Aprop * ve^2; % Corrected thrust formula

% Thrust-to-weight ratio
TWR = FT / FW;

% Display results
fprintf(['\nTHRUST:\n' ...
    'Propeller Diameter:\n  Dprop = 12 x 0.0254\n' ...
    '        = %.4f [m]\n' ...
    'Propeller Radius:\n  rprop = 0.3048/2\n' ...
    '        = %.4f [m]\n' ...
    'Propeller Area:\n  Aprop = pi x rprop\n' ...
    '        = %.3f [m^2]\n' ...
    'Specs Max Thrust:\n  FTmax = 2300 x g / 1000\n' ...
    '        = %.1f [N]\n' ...
    'Exit Velocity:\n  ve = sqrt(2 * FTmax / (rho * A))\n' ...
    '     = %.4f [m/s]\n' ...
    'Thrust:\n  FT = 1/2 * rho * A * ve^2\n' ...
    '     = %.4f [m/s]\n' ...
    'Thrust to Weight Ratio:\n  TWR = FT / FW\n' ...
    '      = %.2f'], ...
    Dprop, rprop, Aprop, FTmax, ve, FT, TWR);

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
fprintf('\nCruise/Stall Ratio: %.2f\n', CSR)

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
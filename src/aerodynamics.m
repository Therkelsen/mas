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

b = 1; % Wing span [m]
c = 0.25; % Chord Length [m]
A = c*b; % Wing Area [m^2]
Vcruise = 25; % Cruise Speed [m/s]
fprintf('\nWing Area = Wing Span x Chord Length\n  A = b x c\n    = %i [m] x %.2f [m]\n    = %.2f [m^2]\n', b, c, A)

Re = (Vcruise * c) / nu;

fprintf(['\nKinematic Viscosity at T = %i °C:\n  ν = %i [m^2/s]' ...
    '\nCruise Speed:\n  Vcruise = %.2f [m/s]\n' ...
    '\nReynolds Number:\n  Re = %i\n'], T, nu, Vcruise, Re)

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
% FL = 1/2*rho*Vcruise^2*A*Cl; % N
FL = FW;
% Vcruise = sqrt((2*FL)/(rho*A*Cl));
Cl = (2*FL)/(rho*Vcruise^2*A);

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
FD = 1/2*rho*Vcruise^2*A*Cd; % Drag Force [N]

fprintf(['\nDRAG\nDrag Force:\n  FD = %.2f [N]\n'], FD)

%% Thrust Stuff
Tmax = 2.3; % Max "thrust" of motor
% Static thrust
Tstatic = Tmax * g; % N, when static. Theoretical number, will be smaller when in flight.
Treq = FD % Required thrust to counteract Drag [N]

fprintf(['\nTHRUST\nMax Thrust Spec:\n  m = %.2f' ...
    '\nStatic Thrust:\n  Tstatic = m x g\n' ...
    '          = %.2f x %.2f\n' ...
    '          = %.2f [N]' ...
    '\nRequired Thrust:\n  Treq = %.2f [N]\n'], ...
    m, m, g, Tstatic, Treq)

% Thrust-to-drag ratio
SRR = Tstatic/Treq;

% Is the thrust enough to counteract drag of drone?
if (Tstatic>Treq)
    fprintf('\nYes, we have enough thrust!')
else
    fprintf('\nNo, we do not have enough thrust!')
end
fprintf('\nT/D Ratio: %.2f\n', SRR)

Preq = Treq * Vcruise; % Required power for level flight [W]

eta_prop = 0.8; % Propellor efficiency
eta_motor = 0.85; % Motor efficiency

Pelec = Preq/(eta_prop * eta_motor); % Electrical power input

fprintf(['\nPOWER\nRequired Power:\n' ...
    '  Preq = Treq x Vcruise\n' ...
    '      = %.2f [N] x %.2f [m/s]\n' ...
    '      = %.2f [W]\n' ...
    'Propeller Efficiency:\n  ηprop = %i%%\n' ...
    'Motor Efficiency:\n  ηmotor = %i%%\n' ...
    'Electrical Power Input:\n' ...
    '  Pelec = Preq/(ηprop * ηmotor)\n' ...
    '        = %.2f [W]' ...
    '\n'], ...
    Treq, Vcruise, Preq, eta_prop*100, eta_motor*100, Pelec)

Vcell = 3.7; % Cell nominal voltage [V]
Vbattery = 4*Vcell; % Battery nominal voltage [V]

fprintf(['\nVOLTAGE\nBattery nominal voltage:\n' ...
    '  Vbattery = 4 x Vcell\n' ...
    '           = 4 x %.2f [V]\n' ...
    '           = %.2f [V]\n'], ...
    Vcell, Vbattery)

Idraw = Pelec/Vbattery; % Required current draw [A]
Imax = 41; % Max current draw of motor [A]

fprintf(['\nCURRENT\nRequired power draw:\n' ...
    '  Idraw = Pelec / Vbattery\n' ...
    '        = %.2f [W] / %.2f [V]\n' ...
    '        = %.2f [A]\n' ...
    'Motor maximum current\n' ...
    '  Imax = %.2f [A]\n'], ...
    Pelec, Vbattery, Idraw, Imax)

DMR = Imax/Idraw;
% Can the motor handle the current?
if (Idraw<Imax)
    fprintf('\nYes, the motor can handle\nthe current draw!')
else
    fprintf('\nNo, the motor can not handle\nthe current draw!')
end
fprintf('\nDraw/Max Ratio: %.2f\n', DMR)

%% Stall Stuff
Clstall = 1.35; % Coefficient of Lift at Stall
Cdstall = 0.05; % Coefficient of Drag at Stall
alphastall = 14; % Degrees

% Stall Speed
Vstall = sqrt((2*m*g)/(rho*A*Clstall));

fprintf(['\nSTALL\nLift Coefficient at Stall:\n  Clstall = %.2f\n' ...
    'Drag Coefficient at Stall:\n  Cldrag = %.2f\n' ...
    'Stall Speed: \n  Vstall = %.2f\n'], Clstall, Cdstall, Vstall)

% Cruise to Stall ratio
CSR = Vcruise/Vstall;

% Is our cruise speed higher than our stall speed?
if (Vcruise>Vstall)
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
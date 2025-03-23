clc; clear; close all;
format compact;

%% Load Required Data
aero = load("aerodynamics.mat");
beam = load("wing_beam_materials.mat");

m = aero.m; % Drone mass [kg]
mcomp = aero.m_comp; % Component mass [kg]
mbeam = beam.mbeam; % Beam mass [kg]
g = aero.g; % Gravity [m/s^2]

mbody = m - mcomp - 2*mbeam; % Available weight left for body [kg]

fprintf(['\nDRONE PARAMETERS\nMass of Drone:\n  m = %.2f [kg]\n' ...
    'Mass of Components:\n  mcomp = %.2f [kg]\n' ...
    'Mass of Beam:\n  mbeam = %g [kg]\n' ...
    'Mass of Body:\n  mbody = m - mcomp - 2 x mbeam\n        = %.2f [kg]\n'], m, mcomp, mbeam, mbody);

% Impact velocity - assuming worst-case free fall from height h
vcruise = aero.vcruise;
h = aero.h; % Fall height [m]
vfall = sqrt(2 * g * h); 
% Horizontal and Vertical turns into a Diagonal velocity
vimpact = sqrt(vcruise^2 + vfall^2); % Impact velocity [m/s]

t = aero.t; % Airfoil Maximum Thickness [m]
Aimpact = t^2; % Impact area [m^2] - Front area of body considering quite a steep wing pitch.

fprintf(['\nIMPACT ANALYSIS\nCruise Speed:\n  vcruise = %i [m/s]\n' ...
    'Fall Height:\n  h = %i [m]\n' ...
    'Gravity:\n  g = %.2f [m/s^2]\n' ...
    'Fall Velocity:\n  vfall = sqrt(2 x g x h)\n        = sqrt(2 x %.2f x %i)\n        = %.2f [m/s]\n' ...
    'Impact Velocity:\n  vimpact = %.2f [m/s]\n' ...
    'Impact Area:\n  Aimpact = %.4f [m^2]\n'], ...
    vcruise, h, g, g, h, vfall, vimpact, Aimpact);

%% Is material light enough?
b = aero.b; % Wingspan [m]
c = aero.c; % Chord Length [m]
A = aero.A; % Wing Area [m^2]
t = aero.t; % Airfoil Maximum Thickness [m]

Vbody = A*t; % Body Volume [m^3]

rhomax = mbody/Vbody; % Maximum density of body material [kg/m^3]

% XPS
rhoxps = 55; % Density of XPS [kg/m^3]

fprintf(['\nMATERIAL DENSITY CHECK\nWingspan:\n  b = %i [m]\n' ...
    'Chord Length:\n  c = %.2f [m]\n' ...
    'Wing Area:\n  A = %.2f [m^2]\n' ...
    'Airfoil Maximum Thickness:\n  t = %.2f [m]\n' ...
    'Body Volume:\n  Vbody = A x t\n        = %.2f [m^3]\n'], ...
    b, c, A, t, Vbody)

fprintf('\nXPS Material Density Check:');
if rhomax > rhoxps
    fprintf('\n  ✓ Low enough density.');
else
    fprintf('\n  ✗ NOT low enough density.');
end

%% Compute Required Stopping Distance and Material Strength
sigma_comp_range = [200E3, 5E6]; % Range of possible material strengths [Pa]
d_stop = (0.5 * m * vimpact^2) ./ (Aimpact * sigma_comp_range);

fprintf(['\nREQUIRED STOPPING DISTANCE\nFor Material Strength %.1f kPa:\n  d_stop = %.4f [m]\n' ...
    'For Material Strength %.1f MPa:\n  d_stop = %.4f [m]\n'], sigma_comp_range(1)/1E3, d_stop(1), sigma_comp_range(2)/1E6, d_stop(2));

%% Check if XPS or Balsa can withstand the impact
Exps = 3.34E9; % Young's Modulus for XPS [Pa]
tens_xps = 250E3; % Tensile Strength [Pa]
comp_xps = 450E3; % Compressive Strength [Pa]

fprintf('\nMATERIAL CHECK\n');
fprintf('XPS Foam: Compressive Strength = %.1f kPa, Required d_stop = %.4f [m]\n', comp_xps/1E3, (0.5 * m * vimpact^2) / (Aimpact * comp_xps));
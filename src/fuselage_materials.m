clc; clear; close all;
format compact;

%% Load Required Data
aero = load("aerodynamics.mat");
beam = load("wing_spar_materials.mat");

m = aero.m; % Drone mass [kg]
mcomp = aero.m_comp; % Component mass [kg]
mbeam = beam.mbeam; % Beam mass [kg]
g = aero.g; % Gravity [m/s^2]

mbody = m - mcomp - 2*mbeam; % Available weight left for body [kg]

fprintf(['\nDRONE PARAMETERS\nMass of Drone:\n  m = %.2f [kg]\n' ...
    'Mass of Components:\n  mcomp = %.2f [kg]\n' ...
    'Mass of Beam:\n  mbeam = %g [kg]\n' ...
    'Mass of Body:\n  mbody = m - mcomp - 2 x mbeam\n        = %.2f [kg]\n'], m, mcomp, mbeam, mbody);

%% Is material light enough?
b = aero.b; % Wingspan [m]
c = aero.c; % Chord Length [m]
A = aero.A; % Wing Area [m^2]
t = aero.t; % Airfoil Maximum Thickness [m]

Vbody = A*t; % Body Volume [m^3]

rhomax = mbody/Vbody % Maximum density of body material [kg/m^3]

% XPS
rhoxps = 55; % Density of XPS [kg/m^3]

fprintf(['\nMATERIAL DENSITY CHECK\nWingspan:\n  b = %i [m]\n' ...
    'Chord Length:\n  c = %.2f [m]\n' ...
    'Wing Area:\n  A = %.2f [m^2]\n' ...
    'Airfoil Maximum Thickness:\n  t = %.2f [m]\n' ...
    'Body Volume:\n  Vbody = A x t\n        = %.2f [m^3]\n' ...
    'Max Body Density:\n  rhomax = mbody/Vbody\n         = %.2f [kg/m^3]\n' ...
    'XPS Density:\n  rhoxps = %i [kg/m^3]\n'], ...
    b, c, A, t, Vbody, rhomax, rhoxps);

fprintf('\nXPS Material Density Check:');
if rhomax > rhoxps
    fprintf('\n  ✓ Low enough density.\n');
else
    fprintf('\n  ✗ NOT low enough density.\n');
end
fprintf(['\nDensity Ratio = rhomax/rhoxps\n' ...
    '              = %.2f\n'], rhomax/rhoxps);
%% Impact velocity - assuming worst-case free fall from height h
vcruise = aero.vcruise;
h = aero.h; % Fall height [m]
vfall = sqrt(2 * g * h); 
vimpact = sqrt(vcruise^2 + vfall^2); % Impact velocity [m/s]
d_stop = 0.05; % Stopping distance on grass field [m]
Aimpact = aero.t^2; % Impact area [m^2] (Assumed impact area)
comp_xps = 450E3; % Compressive Strength of XPS [Pa]

% Calculate the maximum height the material can withstand
h_max = abs((Aimpact * comp_xps * d_stop) / (m * g) - (vcruise^2 / (2 * g)));

% Display result
fprintf('\nMAXIMUM ALTITUDE XPS CAN WITHSTAND\n');
fprintf('Maximum fall height: h = %.2f meters\n', h_max);

%% Save Variables
save("fuselage_materials.mat")
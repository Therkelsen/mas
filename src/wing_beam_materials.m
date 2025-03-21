clc; clear; close all;
format compact;

%% Load Required Data
vars = load("aerodynamics.mat");

%% Step 1: Define Basic Parameters
m = vars.m; % Drone mass [kg]
Vdash = 40; % Dash speed [m/s]
n = 5; % Load Factor (G-force)
g = vars.g; % Gravity [m/s^2]

% Wingspan and Spar Dimensions
b = vars.b; % Wingspan [m]
h = 0.025; % Spar height [m]
d = 0.004; % Spar depth [m]
L = b/2; % Spar length (half wingspan) [m]

fprintf(['\nBASIC PARAMETERS\n' ...
    'Drone Mass:\n  m = %.2f [kg]\n' ...
    'Dash Speed:\n  Vdash = %.2f [m/s]\n' ...
    'Load Factor:\n  n = %.2f [G]\n' ...
    'Gravity:\n  g = %.2f [m/s^2]\n' ...
    'Wingspan:\n  b = %.2f [m]\n' ...
    'Spar Height:\n  h = %.3f [m]\n' ...
    'Spar Depth:\n  d = %.3f [m]\n' ...
    'Spar Length:\n  L = %.2f [m]\n'], m, Vdash, n, g, b, h, d, L);

%% Step 2: Compute Worst-Case Aerodynamic Forces
FW = m * g; % Weight force of drone [N]
FL = n * FW; % Required lift force [N]

p = FL / L; % Lift per unit length of wing [N/m]
c = h / 2; % Distance from Neutral Axis to Extreme Fibers [m]

fprintf(['\nFORCES\nDrone Weight:\n  FW = %.2f [N]\n' ...
    'Required Lift Force:\n  FL = %.2f [N]\n' ...
    'Lift per unit length:\n  p = %.2f [N/m]\n' ...
    'Distance from Neutral Axis:\n  c = %.5f [m]\n'], FW, FL, p, c);

%% Step 3: Calculate Structural Mechanics (Bending Stress)
Mmax = (L^2 * p) / 2; % Maximum bending moment [Nm]
I = 1/12 * L * h^3; % Moment of Inertia [kg·m^2]
sigmamax = (abs(Mmax) * c) / I; % Maximum stress [Pa]

fprintf(['\nSTRESS CALCULATION\nMaximum Moment:\n  M_max = %.4f [Nm]\n' ...
    'Moment of Inertia:\n  I = %.6e [m^4]\n' ...
    'Maximum Stress:\n  σmax = %.4f [Pa]\n       = %.4f [MPa]\n'], Mmax, I, sigmamax, sigmamax/(1E6));

%% Step 4: Material Selection (Balsa Wood)
Ebalsa = 2.1E9; % Young's Modulus for Balsa [Pa]
epsilon_balsa = sigmamax / Ebalsa; % Strain

tens_balsa = 10.5E6; % Tensile Strength [Pa]
comp_balsa = 5E6; % Compressive Strength [Pa]

fprintf(['\nBALSA WOOD\nYoung’s Modulus:\n  Ebalsa = %g [Pa]\nStrain:\n  εbalsa = %.2e\n' ...
    'Tensile Strength:\n  σt = %i [Pa]\n' ...
    'Compressive Strength:\n  σc = %i [Pa]\n'], Ebalsa, epsilon_balsa, tens_balsa, comp_balsa);

% Check if Balsa is strong enough
fprintf('\nBalsa Material Strength Check:');
if tens_balsa > sigmamax
    fprintf('\n  ✓ Strong enough in tension.');
else
    fprintf('\n  ✗ NOT strong enough in tension.');
end
if comp_balsa > sigmamax
    fprintf('\n  ✓ Strong enough in compression.');
else
    fprintf('\n  ✗ NOT strong enough in compression.');
end
fprintf('\nStrength Ratios:\n  Tensile = %.2f\n  Compressive = %.2f\n', tens_balsa/sigmamax, comp_balsa/sigmamax);

%% Step 5: Material Selection (XPS Foam)
Exps = 3.34E9; % Young's Modulus for XPS [Pa]
epsilon_xps = sigmamax / Exps;

tens_xps = 250E3; % Tensile Strength [Pa]
comp_xps = 450E3; % Compressive Strength [Pa]

fprintf(['\nEXTRUDED POLYSTYRENE FOAM\nYoung’s Modulus:\n  Exps = %g [Pa]\nStrain:\n  εxps = %.2e\n' ...
    'Tensile Strength:\n  σt = %i [Pa]\n' ...
    'Compressive Strength:\n  σc = %i [Pa]\n'], Exps, epsilon_xps, tens_xps, comp_xps);
% Check if XPS is strong enough
fprintf('\nXPS Material Strength Check:');
if tens_xps > sigmamax
    fprintf('\n  ✓ Strong enough in tension.');
else
    fprintf('\n  ✗ NOT strong enough in tension.');
end
if comp_xps > sigmamax
    fprintf('\n  ✓ Strong enough in compression.');
else
    fprintf('\n  ✗ NOT strong enough in compression.');
end
fprintf('\nStrength Ratios:\n  Tensile = %.2f\n  Compressive = %.2f\n', tens_xps/sigmamax, comp_xps/sigmamax);

%% Step 6: Export Data
save("wing_beam_materials.mat")
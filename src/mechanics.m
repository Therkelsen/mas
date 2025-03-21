clc; clear; close all;
format compact;

%% Mechanics
% load aerodynamics.mat
% requires you to have run aerodynamics.m
vars = load("aerodynamics.mat");

m = vars.m; % Drone mass [kg]
Vdash = 40; % Dash speed [m/s]
n = 5; % Load Factor, i.e. G-forces
g = vars.g; % Gravity

% Force on the Wing Spars
FW = m*g; % Weight force of drone
FL = n*FW; % Required lift force

fprintf(['\nFORCES\nDrone Weight:\n' ...
    '  FW = m * g\n     = %.2f [kg] * %.2f [m/s^2]\n     = %.2f [N]\n' ...
    'Required Lift Force:\n' ...
    '  FL = n * FW\n     = %.2f * %.2f [N]\n     = %.2f [N]\n'], ...
    m, g, FW, n, FW, FL);

% Wingspan
b = vars.b; % m
% Spar height
h = 0.025; % m
% Spar depth
d = 0.004; % m
% Spar length
L = b/2; % m
% Cross-sectional area of airfoil (assuming rectangle)
A = L*h;

% Stress = Tensile Strength
sigma = FL/A; % Pa (N/m²)

%% Cantilever Stuff
% Line Pressure Load aka
% Lift per unit length of wing
p = FL/L; % N/m (Force distributed along the length)

% Distance from Neutral Axis
% to Extreme Fibers
c = h/2; % m

% Moment of Inerta
I = 1/12*L*h^3; % kg/m^2

fprintf(['\nDIMENSIONS\nSpar heigh:\n  h = %.2f [m]\n' ...
    'Spar Length:\n  L = %.2f [m]\n' ...
    'Lift per Unit Length of Wing:\n' ...
    '  p = FL / L\n    = %.2f [N/m]\n' ...
    'Distance from Neutral Axis to Extreme Fibers:\n  c = %.2f [m]\n' ...
    'Moment of Inertia:\n' ...
    '  I = 1/12 * L * h^3 = %i [kg·m^2]\n'], ...
    h, L, p, c, I);

% Maximum Moment
Mmax = (L^2*p)/2; % Nm

% Maximum Stress
sigmamax = (abs(Mmax)*c)/I; % Pa

fprintf(['\nSTRESS CALCULATION\nMaximum Moment (Mmax):\n' ...
    '  M_max = (L^2 * p) / 2 = %.4f [Nm]\n' ...
    'Maximum Stress (σmax):\n' ...
    '  σmax = (abs(Mmax) * c) / Iedgezz\n' ...
    '       = %.4f [Pa]\n' ...
    '       = %.4f [MPa]\n'], ...
    Mmax, sigmamax, sigmamax/(1E6));

%% Balsa
% Young's Modulus for Balsa
Ebalsa = 2.1E9;

% Strain
epsilon_balsa = sigma/Ebalsa;

% Tensile Strength
tens_balsa = 10.5E6;
% Compressive Strength
comp_balsa = 5E6;

fprintf(['\nBALSA WOOD\nYoungs Modulus:\n  Ebalsa = %i\n' ...
    'Balsa Strain:\n  εbalsa = %i\n' ...
    'Tensile Strength = %i [MPa]\n' ...
    'Compressive Strength = %i [MPa]\n'], Ebalsa, epsilon_balsa, tens_balsa, comp_balsa)

% Will the material break under stress?
% Max Stress vs Material Tensile Strength
if (tens_balsa>sigmamax)
    fprintf('\nYes, our material is strong enough to withstand tensile stress.')
else
    fprintf('\nNo, our is not strong enough to withstand tensile stress.')
end
fprintf('\nMax Stress vs Material Tensile Strength Ratio: %.2f\n', tens_balsa/sigmamax)

% Max Stress vs Material Compressive Strength
if (comp_balsa>sigmamax)
    fprintf('\nYes, our material is strong enough to withstand compressive stress.')
else
    fprintf('\nNo, our is not strong enough to withstand compressive stress.')
end
fprintf('\nMax Stress vs Material Compressive Strength Ratio: %.2f\n', comp_balsa/sigmamax)

%% XPS (Extruded Polystyrene Foam) 
% Young's Modulus for XPS
Exps = 3.34E9;
% Strain
epsilon_xps = sigma/Exps;

%Compressive strength of 250 kPa, Tensile strength of 450 kPa

% Tensile Strength
tens_xps = 250E3;
% Compressive Strength
comp_xps = 450E3;

fprintf(['\nEXTRUDED POLYSTYRENE FOAM\nYoungs Modulus:\n  Exps = %i\n' ...
    'XPS Strain:\n  εxps= %i\n' ...
    'Tensile Strength = %i [MPa]\n' ...
    'Compressive Strength = %i [MPa]\n'], Exps, epsilon_xps, tens_xps, comp_xps)

% Will the material break under stress?
% Max Stress vs Material Tensile Strength
if (tens_xps>sigmamax)
    fprintf('\nYes, our material is strong enough to withstand tensile stress.')
else
    fprintf('\nNo, our is not strong enough to withstand tensile stress.')
end
fprintf('\nMax Stress vs Material Tensile Strength Ratio: %.2f\n', tens_xps/sigmamax)

% Max Stress vs Material Compressive Strength
if (comp_xps>sigmamax)
    fprintf('\nYes, our material is strong enough to withstand compressive stress.')
else
    fprintf('\nNo, our is not strong enough to withstand compressive stress.')
end
fprintf('\nMax Stress vs Material Compressive Strength Ratio: %.2f\n', comp_xps/sigmamax)

%% Export Variables
save("mechanics.mat")
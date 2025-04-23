clc; clear; close all;
format compact;

%% Load Required Data
aero = load("aerodynamics.mat");
spar = load("wing_spar_materials.mat");

m = aero.m; % Drone mass [kg]
m_comp = aero.m_comp; % Component mass [kg]
m_spar = spar.mspar; % Beam mass [kg]
g = aero.g; % Gravity [m/s^2]

m_body = m - m_comp - m_spar; % Available weight left for body [kg]

fprintf('\nDRONE MASS BUDGET\n');
fprintf('  Total Drone Mass:\n    m = %.2f [kg]\n', m);
fprintf('  Component Mass:\n    m_comp = %.2f [kg]\n', m_comp);
fprintf('  Beam Mass:\n    m_spar = %.2f [kg]\n', m_spar);
fprintf('  Remaining for Body:\n    m_body = %.2f [kg]\n', m_body);

%% Is material light enough?
h_body = 0.08; % [m] outer diameter
l_body = 0.5;  % [m] length of body
t_body = 0.005; % [m] wall thickness

r_outer = h_body / 2;
r_inner = r_outer - t_body;

V_fuselage = pi * l_body * r_outer^2;
V_fuselage_shell = pi * l_body * (r_outer^2 - r_inner^2);
rho_max = m_body/V_fuselage; % Max allowable average body density [kg/m^3]

% XPS
rho_XPS = 55; % Density of XPS [kg/m^3]

fprintf('\nMATERIAL DENSITY CHECK\n');
fprintf('  Fuselage Volume:\n    V_fuselage = %.6f [m^3]\n', V_fuselage);
fprintf('  Max Allowable Density:\n    rho_max = %.2f [kg/m^3]\n', rho_max);

fprintf('\nXPS Density:\n   rho_XPS    = %.2f [kg/m^3]\n', rho_XPS);
fprintf('Density Ratio (XPS) = %.2f\n', rho_max / rho_XPS);
if rho_XPS < rho_max
    fprintf('   ✓ XPS meets density requirement.\n');
else
    fprintf('   ✗ XPS too heavy.\n');
end

m_body_XPS = rho_XPS * V_fuselage; % XPS body mass [kg]
fprintf('  XPS Body Volume:\n    V_fuselage_shell = %.6f [m^3]\n', V_fuselage);
fprintf('  XPS Body Mass:\n    m_body_XPS = %.2f [kg]\n', m_body_XPS);

%% Add PLA Shell Material Check
fprintf('\n\nPLA SHELL CHECK\n');

rho_max = m_body/V_fuselage_shell; % Max allowable average body density [kg/m^3]
rho_PLA = 1240; % Density of PLA [kg/m^3]

fprintf('\nPLA Density:\n   rho_PLA    = %.2f [kg/m^3]\n', rho_PLA);
fprintf('Density Ratio (PLA) = %.2f\n', rho_max / rho_PLA);
if rho_PLA < rho_max
    fprintf('  ✓ PLA meets density requirement.\n');
else
    fprintf('  ✗ PLA too heavy.\n');
end

m_body_PLA = rho_PLA * V_fuselage_shell; % PLA shell mass [kg]
fprintf('  PLA Shell Volume:\n    V_fuselage = %.6f [m^3]\n', V_fuselage_shell);
fprintf('  PLA Shell Mass:\n    m_body_PLA = %.2f [kg]\n', m_body_PLA);

%% Save Variables
save("fuselage_materials.mat")

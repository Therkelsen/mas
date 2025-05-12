clc; clear; close all;
format compact;

%% Load Required Data
aero = load("aerodynamics.mat");
spar = load("wing_spar_materials.mat");

m = aero.m % Drone mass [kg]
m_comp = aero.m_comp % Component mass [kg]
m_spar = spar.mspar % Beam mass [kg]
g = aero.g; % Gravity [m/s^2]

m_body = m - m_comp - m_spar; % Available weight left for body [kg]

fprintf('\nDRONE MASS BUDGET\n');
fprintf('  Total Drone Mass:\n    m = %.2f [kg]\n', m);
fprintf('  Component Mass:\n    m_comp = %.2f [kg]\n', m_comp);
fprintf('  Beam Mass:\n    m_spar = %.2f [kg]\n', m_spar);
fprintf('  Remaining for Body:\n    m_body = %.2f [kg]\n', m_body);

%% Is material light enough?
h_body = 0.08; % [m] outer diameter
l_body = 0.45; % [m] length of body
t_body = 0.005; % [m] wall thickness

r_outer = h_body / 2;
r_inner = r_outer - t_body;

V_fuse = pi * l_body * (r_outer^2 - r_inner^2);
rho_max = 0.79 / V_fuse; % Max allowable average shell material density [kg/m^3]

fprintf('\nMATERIAL DENSITY CHECK\n');
fprintf('  Fuselage Shell Volume:\n    V_fuse = %.6f [m^3]\n', V_fuse);
fprintf('  Max Allowable Shell Density:\n    rho_max = %.2f [kg/m^3]\n', rho_max);

%% XPS Shell
rho_XPS = 55; % Density of XPS [kg/m^3]

fprintf('\nXPS SHELL CHECK\n');
fprintf('XPS Density:\n    rho_XPS = %.2f [kg/m^3]\n', rho_XPS);
fprintf('Density Ratio (XPS) = %.2f\n', rho_max / rho_XPS);
if rho_XPS < rho_max
    fprintf('   ✓ XPS meets shell density requirement.\n');
else
    fprintf('   ✗ XPS too heavy.\n');
end

m_XPS = rho_XPS * V_fuse; % XPS shell mass [kg]
fprintf('XPS Shell Mass:\n    m_XPS = %.2f [kg]\n', m_XPS);

%% PLA Shell
rho_PLA = 1240; % Density of PLA [kg/m^3]

fprintf('\nPLA SHELL CHECK\n');
fprintf('PLA Density:\n    rho_PLA = %.2f [kg/m^3]\n', rho_PLA);
fprintf('Density Ratio (PLA) = %.2f\n', rho_max / rho_PLA);
if rho_PLA < rho_max
    fprintf('   ✓ PLA meets shell density requirement.\n');
else
    fprintf('   ✗ PLA too heavy.\n');
end

% m_body_PLA = rho_PLA * V_fuselage_shell; % PLA shell mass [kg]
m_body_PLA = 0.159;
fprintf('  PLA Shell Volume:\n    V_fuselage = %.6f [m^3]\n', V_fuselage_shell);
fprintf('  PLA Shell Mass:\n    m_body_PLA = %.2f [kg]\n', m_body_PLA);

%% Save Variables
save("fuse_materials.mat")

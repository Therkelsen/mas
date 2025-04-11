clc; clear; close all;
format compact;

%% Load Required Data
aero = load("aerodynamics.mat");
beam = load("wing_spar_materials.mat");
body = load("fuselage_materials.mat");

m = aero.m; % Drone mass [kg]
n = beam.n; % Load Factor (G-force)
g = aero.g; % Gravity [m/s^2]

FT = aero.FT; % Motor Thrust (12x6 prop) [N]

%% Pulling Force
Fpull = n*g*m;

fprintf(['\nFORCE PULLING ON MOUNT\n' ...
    'Force in vertical maneuver:\n' ...
    'Fpull  = n x g x m\n       = %.2f [N]\n'], ...
    Fpull);

%% Torque
% Propeller specifications (12x6 inch prop)
Dstator = 0.02; % Stator Diameter [m]

% Total distance from motor center to prop hub
r = Dstator / 2; % in meters

% Calculate torque due to motor thrust force (Torque = Force * Distance)
tau = FT * r;

fprintf(['\nTORQUE ON MOUNT\n' ...
    'Torque at full throttle:\n' ...
    '  tau = FT x Dstator/2 =\n' ...
    '      = %g [Nm]\n'], ...
    tau);

%% Material Properties (Rubber, CFRP, GFRP)
% Rubber (EPDM)
rubber_tensile_strength = 20e6; % Tensile Strength [Pa] = 20 MPa
rubber_compressive_strength = 15e6; % Compressive Strength [Pa] = 15 MPa
rubber_shear_strength = 15e6; % Shear Strength [Pa] = 15 MPa

% CFRP (Carbon Fiber Reinforced Polymer)
CFRP_tensile_strength = 1900e6; % Tensile Strength [Pa] = 1200 MPa
CFRP_compressive_strength = 1300e6; % Compressive Strength [Pa] = 700 MPa
CFRP_shear_strength = 80e6; % Shear Strength [Pa] = 80 MPa

% GFRP (Glass Fiber Reinforced Polymer)
GFRP_tensile_strength = 800e6; % Tensile Strength [Pa] = 800 MPa
GFRP_compressive_strength = 400e6; % Compressive Strength [Pa] = 400 MPa
GFRP_shear_strength = 50e6; % Shear Strength [Pa] = 50 MPa

%% Calculate Stress on Mount

% Assume the torque is applied over a cylindrical
% area on the motor mount (simplification)
mount_radius = r; % Example mount radius [m]
% (this can be adjusted based on the design)

% Area of the cylindrical motor mount
mount_area = pi * (mount_radius^2);

% Calculate stress due to torque (Stress = Torque / (Radius * Area))
stress_torque = tau / (mount_radius * mount_area);

% Calculate stress due to pulling force (Shear Stress = Force / Area)
stress_force = Fpull / mount_area;

fprintf(['\nMATERIAL STRESS CHECK\n' ...
    'Rubber (EPDM)\n' ...
    '  Tensile strength:\n' ...
    '    sigma_t = %g [Pa]\n' ...
    '  Compressive strength:\n' ...
    '    sigma_c = %g [Pa]\n' ...
    '  Shear strength:\n' ...
    '    sigma_s = %g [Pa]\n' ...
    'CFRP\n' ...
    '  Tensile strength:\n' ...
    '    sigma_t = %g [Pa]\n' ...
    '  Compressive strength:\n' ...
    '    sigma_c = %g [Pa]\n' ...
    '  Shear strength:\n' ...
    '    sigma_s = %g [Pa]\n' ...
    'GFRP\n' ...
    '  Tensile strength:\n' ...
    '    sigma_t = %g [Pa]\n' ...
    '  Compressive strength:\n' ...
    '    sigma_c = %g [Pa]\n' ...
    '  Shear strength:\n' ...
    '    sigma_s = %g [Pa]\n'], ...
    rubber_tensile_strength, rubber_compressive_strength, rubber_shear_strength, ...
    CFRP_tensile_strength, CFRP_compressive_strength, CFRP_shear_strength, ...
    GFRP_tensile_strength, GFRP_compressive_strength, GFRP_shear_strength);

%% Print Material Stress Check

fprintf(['\nMATERIAL STRESS CHECK\n' ...
    '  Stress from Torque = %.2f [Pa]\n' ...
    '  Stress from Force = %.2f [Pa]\n'], ...
    stress_torque, stress_force);

%% Material Strength Comparison and Pass/Fail Check

fprintf('\nMATERIAL STRENGTH CHECK:\n');

% Rubber (EPDM) strength checks
fprintf('  Rubber (EPDM):\n');
if stress_torque > rubber_shear_strength
    fprintf('    ✗ Cannot handle the shear stress.\n');
else
    fprintf('    ✓ Can handle the shear stress.\n');
end
fprintf('       Shear Strength Ratio: %.2f\n', rubber_shear_strength / stress_torque);

if stress_force > rubber_tensile_strength
    fprintf('    ✗ Cannot handle the tensile stress.\n');
else
    fprintf('    ✓ Can handle the tensile stress.\n');
end
fprintf('       Tensile Strength Ratio: %.2f\n', rubber_tensile_strength / stress_force);

if stress_force > rubber_compressive_strength
    fprintf('    ✗ Cannot handle the compressive stress.\n');
else
    fprintf('    ✓ Can handle the compressive stress.\n');
end
fprintf('       Compressive Strength Ratio: %.2f\n', rubber_compressive_strength / stress_force);

% CFRP strength checks
fprintf('  CFRP:\n');
if stress_torque > CFRP_shear_strength
    fprintf('    ✗ Cannot handle the shear stress.\n');
else
    fprintf('    ✓ Can handle the shear stress.\n');
end
fprintf('       Shear Strength Ratio: %.2f\n', CFRP_shear_strength / stress_torque);

if stress_force > CFRP_tensile_strength
    fprintf('    ✗ Cannot handle the tensile stress.\n');
else
    fprintf('    ✓ Can handle the tensile stress.\n');
end
fprintf('       Tensile Strength Ratio: %.2f\n', CFRP_tensile_strength / stress_force);

if stress_force > CFRP_compressive_strength
    fprintf('    ✗ Cannot handle the compressive stress.\n');
else
    fprintf('    ✓ Can handle the compressive stress.\n');
end
fprintf('       Compressive Strength Ratio: %.2f\n', CFRP_compressive_strength / stress_force);

% GFRP strength checks
fprintf('  GFRP:\n');
if stress_torque > GFRP_shear_strength
    fprintf('    ✗ Cannot handle the shear stress.\n');
else
    fprintf('    ✓ Can handle the shear stress.\n');
end
fprintf('       Shear Strength Ratio: %.2f\n', GFRP_shear_strength / stress_torque);

if stress_force > GFRP_tensile_strength
    fprintf('    ✗ Cannot handle the tensile stress.\n');
else
    fprintf('    ✓ Can handle the tensile stress.\n');
end
fprintf('       Tensile Strength Ratio: %.2f\n', GFRP_tensile_strength / stress_force);

if stress_force > GFRP_compressive_strength
    fprintf('    ✗ Cannot handle the compressive stress.\n');
else
    fprintf('    ✓ Can handle the compressive stress.\n');
end
fprintf('       Compressive Strength Ratio: %.2f\n', GFRP_compressive_strength / stress_force);

save('motor_mount_materials.mat')
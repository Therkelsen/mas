clc; clear; close all;
format compact;

aero = load('aerodynamics.mat');
fuse = load('fuselage_materials.mat');
moto = load('motor_mount_materials.mat');
wing = load('wing_spar_materials.mat');

%% Stability
swp_ang = 20; % sweep angle [deg]
swp_ang_rad = deg2rad(swp_ang); % [rad]

xnp = 0.05 + 1/4 * aero.c * (1 - tan(swp_ang_rad)); % Neutral Point [m]

% Masses and positions (example)
masses = [aero.m_motor, aero.m_payload, ...
          aero.m_battery, aero.m_battery];   % Masses of objects [kg]

x_motor = 0;
x_body = 0;
x_payload = 0.045;
x_battery = xnp;

l_motor = 0.045; % Motor length [m]
l_body = 0.5; % Body length [m]
l_payload = 0.09; % Payload length [m]
l_battery = 0.055; % Battery length [m]

positions = [x_motor+l_motor/2, x_payload+l_payload/2, ...
             x_battery+l_battery/2, x_body+l_battery/2]; % Positions of objects [m]

% Calculate the center of mass
total_mass = sum(masses);            % Total mass
weighted_sum = sum(masses .* positions); % Weighted sum of masses and positions
xcg = weighted_sum / total_mass;     % Center of mass
sm = (xnp - xcg)/aero.c;             % Stability Margin (between 0.05 and 0.4 is acceptable)

lh = l_body - xnp; % horizontal tail moment arm [m]
ch = 0.125;        % chord length of horizontal rudder [m]
bh = 0.25;         % wingspan of horizontal rudder [m]
Ah = ch*bh;        % Area of horizontal rudder [m^2]

Vh = (Ah*lh)/(aero.Awing*aero.c); % Horizontal Tail Volume Coeff

lv = lh;            % vertical tail moment arm [m]
cv = ch*0.75;       % chord length of vertical rudder [m]
bv = ch*0.75;       % wingspan of vertical rudder [m]
Av = cv*bv;         % Area of vertical rudder [m^2]

Vv = (Av*lv)/(aero.Awing*aero.b); % vertical Tail Volume Coeff

gamma = 3;          % Dihedral angle

B = lv/aero.b * gamma/aero.Cl; % Blaine Stability (<5 unstable, =5 neutral, >5 stable)

% Output results using fprintf
fprintf('\nSTABILITY PARAMETERS:\n');
fprintf('Neutral Point:\n  xnp = %.4f [m]\n', xnp);
fprintf('Center of Mass:\n  xcg = %.4f [m]\n', xcg);
fprintf('Stability Margin:\n  sm = %.4f\n', sm);
fprintf('Horizontal Tail Moment Arm:\n  lh = %.4f [m]\n', lh);
fprintf('Chord Length of Horizontal Rudder:\n  ch = %.4f [m]\n', ch);
fprintf('Wingspan of Horizontal Rudder:\n  bh = %.4f [m]\n', bh);
fprintf('Area of Horizontal Rudder:\n  Ah = %.4f [m^2]\n', Ah);
fprintf('Horizontal Tail Volume Coefficient:\n  Vh = %.4f\n', Vh);

fprintf('\nVertical Tail Parameters:\n');
fprintf('Vertical Tail Moment Arm:\n  lv = %.4f [m]\n', lv);
fprintf('Chord Length of Vertical Rudder:\n  cv = %.4f [m]\n', cv);
fprintf('Wingspan of Vertical Rudder:\n  bv = %.4f [m]\n', bv);
fprintf('Area of Vertical Rudder:\n  Av = %.4f [m^2]\n', Av);
fprintf('Vertical Tail Volume Coefficient:\n  Vv = %.4f\n', Vv);

fprintf('\nBlaine Stability:\n  B = %.4f\n', B);

save('stability.mat');

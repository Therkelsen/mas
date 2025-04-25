clc; clear; close all;
format compact;

aero = load('aerodynamics.mat');
wing = load('wing_spar_materials.mat');
fuse = load('fuselage_materials.mat');

%% Stability
fprintf('\nPITCH STABILITY:\n');

Lambda_LE = 20; % sweep angle [deg]
fprintf('  Wing Sweep Angle:\n    Lambda_LE = %i°\n', Lambda_LE)

Lambda_LE_rad = deg2rad(Lambda_LE); % [rad]

x_airfoil = 0.09;
x_np = x_airfoil + 1/4 * aero.c * (1 - tan(Lambda_LE_rad)); % Neutral Point [m]

l_body = fuse.l_body; % Body length [m]
lh = l_body - x_np;   % horizontal tail moment arm [m]
ch = 0.125;           % chord length of horizontal rudder [m]
bh = 0.35;            % wingspan of horizontal rudder [m]
Ah = ch*bh;           % Area of horizontal rudder [m^2]

V_horizontal = Ah*aero.h_spar;
m_horizontal = V_horizontal * wing.rhobalsa;

Vh = (Ah*lh)/(aero.Awing*aero.c); % Horizontal Tail Volume Coeff

fprintf('\nHORIZONTAL RUDDER\n')
fprintf('  Horizontal Tail Moment Arm:\n    lh = %.4f [m]\n', lh);
fprintf('  Chord Length of Horizontal Rudder:\n    ch = %.4f [m]\n', ch);
fprintf('  Wingspan of Horizontal Rudder:\n    bh = %.4f [m]\n', bh);
fprintf('  Area of Horizontal Rudder:\n    Ah = %.4f [m^2]\n', Ah);
fprintf('  Horizontal Tail Volume Coefficient:\n    Vh = %.4f -> ', Vh);

if Vh < 0.3
    disp('Too Small');
elseif Vh >= 0.3 && Vh <= 0.6
    disp('Within Acceptable Range');
elseif Vh > 0.6
    disp('Too Large');
else
    disp('Unknown stability margin');
end

lv = lh;            % vertical tail moment arm [m]
cv = ch*0.55;       % chord length of vertical rudder [m]
bv = ch*0.85;       % Wingspan of vertical rudder [m]
Av = cv*bv;         % Area of vertical rudder [m^2]

V_vertical = Av*aero.h_spar;
m_vertical = V_vertical * wing.rhobalsa;

Vv = (Av*lv)/(aero.Awing*aero.b); % vertical Tail Volume Coeff

fprintf('\nVERTICAL RUDDER\n')
fprintf('  Vertical Tail Moment Arm:\n    lv = %.4f [m]\n', lv);
fprintf('  Chord Length of Vertical Rudder:\n    cv = %.4f [m]\n', cv);
fprintf('  Wingspan of Vertical Rudder:\n    bv = %.4f [m]\n', bv);
fprintf('  Area of Vertical Rudder:\n    Av = %.4f [m^2]\n', Av);
fprintf('  Vertical Tail Volume Coefficient:\n    Vv = %.4f -> ', Vv);

if Vv < 0.02
    disp('Too Small');
elseif Vv >= 0.02 && Vv <= 0.05
    disp('Within Acceptable Range');
elseif Vv > 0.05
    disp('Too Large');
else
    disp('Unknown stability margin');
end

fprintf('  Neutral Point:\n    x_np = %.4f [m]\n', x_np);

m_rudder = m_vertical*m_horizontal;

m_airfoil = aero.Vwing * fuse.rho_XPS - wing.mspar;

% Masses of objects [kg]
M_vec = [fuse.m_body_PLA, aero.m_motor, ...
          aero.m_payload, m_airfoil, ...
          aero.m_battery, m_rudder];

l_motor = 0.045; % Motor length [m]
l_payload = 0.06; % Payload length [m]
l_battery = 0.055; % Battery length [m]
l_rudder = ch; % Rudder Chord Length

% Lengths of objects
l_vec = [l_body, l_motor, l_payload, aero.c * cos(Lambda_LE_rad), l_battery, l_rudder];

h_body = fuse.h_body; % Body height 80 mm
h_motor = 0.035; % Motor height 35 mm
h_payload = 0.06; % Payload height 6 mm
h_battery = 0.0075; % Battery height 7.5 mm
h_rudder = bv; % Rudder height
h_vec = [h_body, h_motor, h_payload, aero.h_spar, h_battery, h_rudder];

y_body = 0; % Body starting position
y_motor = h_body/2 - h_motor/2; % Motor starting position
y_payload = h_battery; % Payload starting position
y_battery = 0; % Battery starting position
y_airfoil = h_body; % Airfoil starting position
y_rudder = h_body/2; % Rudder starting position

y_vec = [y_body, y_motor, y_payload, y_airfoil, y_battery, y_rudder];

% Positions of objects [m]
x_body = 0; % Body x-position [m]
x_motor = 0; % Motor x-position [m]
x_payload = 0.045; % Payload x-position [m]
x_battery = 0; % Battery x-position [m]
x_rudder = l_body + l_rudder/2; % Wing position

x_vec = [x_body + l_body/2, ...
         x_motor + l_motor/2, ...
         x_payload + l_payload/2, ...
         x_airfoil + aero.c/2, ...
         x_battery + l_battery/2, ...
         x_rudder];

mass_names = {'Body', 'Motor', 'Payload', 'Airfoil', 'Battery', 'Wing'};

% Calculate the center of mass (x_CoM = 1/M * sum(m_i * x_i)
x_cg = sum(M_vec .* x_vec) / sum(M_vec); % Center of mass [m]
fprintf('  Center of Mass:\n    x_cg = %.4f [m]\n', x_cg);

plot_mass_stack(M_vec, x_vec, y_vec, l_vec, h_vec, mass_names, x_np, x_cg);

SM = (x_np - x_cg)/aero.c;             % Stability Margin (between 0.05 and 0.4 is acceptable)
fprintf('  Pitch Stability Margin:\n    SM = %.4f -> ', SM);


if SM <= -0.4
    disp('Strongly Unstable');
elseif SM > -0.4 && SM <= -0.05
    disp('Weakly Unstable');
elseif SM > -0.05 && SM < 0.05
    disp('Neutral');
elseif SM >= 0.05 && SM < 0.4
    disp('Weakly Stable');
elseif SM >= 0.4
    disp('Strongly Stable');
else
    disp('Unknown stability margin');
end

% Dihedral Sizing Criteria – Spiral Stability
gamma = 3; % Dihedral angle

B = lv/aero.b * gamma/aero.Cl; % Blaine Stability (<5 unstable, =5 neutral, >5 stable)

fprintf('\nBLAINE STABILITY:\n');
fprintf('  Dihedral (Upward) Angle:\n    gamma = %i°\n', gamma)
fprintf('  Blaine Stability:\n    B = %.4f -> ', B);

if B < 5
    disp('Spirally Unstable');
elseif B == 5
    disp('Spirally Neutral');
elseif B > 5
    disp('Spirally Stable');
else
    disp('Unknown Blaine Stability');
end

% Dihedral Sizing – Roll Control
VvB = Vv*B;
fprintf('  Roll Control:\n    VvB = %.4f -> ', VvB)

if VvB < 0.1
    disp('Too Low');
elseif VvB >= 0.1 && VvB <= 0.2
    disp('Within Acceptable Range');
elseif VvB > 0.2
    disp('Too High');
else
    disp('Unknown Roll Controllability');
end

%% Save Params
save('stability.mat');

%% Visualization
function plot_mass_stack(M_vec, x_vec, y_vec, l_vec, h_vec, mass_names, x_np, x_cg)
    figure;
    hold on;
    colors = lines(length(M_vec)); % Create distinct colors for each mass
    
    y_pos = 0; % Start at the bottom of the plot for the first mass
    h = zeros(1, length(M_vec) + 2); % Array to store handles to plot objects for the legend
    
    for i = 1:length(M_vec)
        % Calculate the position on the x-axis for each mass (center - half length)
        x_start = x_vec(i) - l_vec(i)/2; % Starting position for each mass on x-axis 
        
        % Plot each mass as a rectangle (same height but stacked on y-axis)
        rectangle('Position', [x_start, y_vec(i), l_vec(i), h_vec(i)], ...
                  'FaceColor', colors(i,:), 'EdgeColor', 'k', 'LineWidth', 1.5);
        
        % Create invisible plot handles for the legend for each mass
        h(i) = plot(NaN, NaN, 's', 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 10);
        
        % Increment y position for the next mass (add the height of the current mass)
        y_pos = y_pos + h_vec(i); 
    end
    
    % Format the mass names to include the mass value next to the name
    mass_names_with_values = cell(1, length(M_vec));
    for i = 1:length(M_vec)
        mass_names_with_values{i} = sprintf('%s (%.3f kg)', mass_names{i}, M_vec(i));
    end
    
    % xlim([-0.05, 0.55])
    % ylim([-0.05, 0.2])

    % Calculate the earliest and latest object positions
    x_min = min(x_vec + min(l_vec)) - 0.1;  % 0.1 before the earliest object
    x_max = max(x_vec + min(l_vec)/16) + 0.1;  % 0.1 after the last object
    
    % Apply the same limits to both axes
    xlim([x_min, x_max]);
    ylim([x_min, x_max]);

    % Plot vertical lines for x_np and x_cg
    h(length(M_vec) + 1) = line([x_np, x_np], ylim, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 2); % Red line for x_np
    h(length(M_vec) + 2) = line([x_cg, x_cg], ylim, 'Color', 'b', 'LineStyle', '--', 'LineWidth', 2); % Blue line for x_cg
    
    % Labeling the plot
    xlabel('X-Position [m]');
    ylabel('Y-Position [m]');
    title('Stacked Mass Distribution');
    grid on;
    
    % Create a legend with mass names and their corresponding values, as well as x_np and x_cg
    legend([h(1:length(M_vec)), h(length(M_vec)+1), h(length(M_vec)+2)], ...
           [mass_names_with_values, {'x_{np}', 'x_{cg}'}], 'Location', 'best');
    
    hold off;
end
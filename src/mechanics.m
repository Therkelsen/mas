clc; clear; close all;
format compact;

% Force on the Wing Spars
FL = 50.12; % N
FW = 14.73; % N

% % Net Load on both wings
F_wing = FL - FW % N
% % Per wing
F_spar = F_wing/2

% Wing Span
b = 1.30; % m
% Airfoil Maximum Thickness
h = 0.03; % m
% Cross-sectional area of airfoils
A = b*h

% Stress
sigma = F_wing/A % Pa

% Young's Modulus for Plywood
E = 10E9;

% Strain
epsilon = sigma/E

% Length of beam
L = b/2 % m
% Line Pressure Load aka
% Lift per unit length of wing
p = FL/b % N/m

% Distance from Neutral Axis
% to Extreme Fibers
c = h/2 % m

% Moment of Inerta
I = 1/12*b*h^3 % kg/m^2

% Maximum Moment
M_max = (L^2*p)/2 % Nm

% Maximum Stress
sigma_max = (abs(M_max)*c)/I % Pa
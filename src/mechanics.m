clc; clear; close all;
format compact;

% Force on the Wing Spars
L = 50.12; % N
W = 14.73; % N

% % Net Load on both wings
F_wing = L - W % N
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
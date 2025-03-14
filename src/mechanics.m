clc; clear; close all;
format compact;

% load aerodynamics.mat
% requires you to have run aerodynamics.m
vars = load("aerodynamics.mat");

% Force on the Wing Spars
FL = vars.FL % N

% Wing Span
b = vars.b % m
% Airfoil Maximum Thickness
h = 0.025; % m
% Length of beam
L = b/2 % m
% Cross-sectional area of airfoil (assuming rectangle)
A = L*h

% Stress
sigma = FL/A % Pa (N/m²)

% Young's Modulus for Plywood
E = 10E9;

% Strain
epsilon = sigma/E

% Line Pressure Load aka
% Lift per unit length of wing
p = FL/L % N/m (Force distributed along the length)

% Distance from Neutral Axis
% to Extreme Fibers
c = h/2 % m

% Moment of Inerta
I = 1/12*L*h^3 % kg/m^2

% Maximum Moment
M_max = (L^2*p)/2 % Nm

% Maximum Stress
sigma_max = (abs(M_max)*c)/I % Pa

% Export Variables
save("mechanics.mat")
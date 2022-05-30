clc 
clear 
close all 

% Reference natural frequency 
omega_ref = 2*pi*0.4;

% Pinned-Pinned 
P1 = fminbnd(@(P) optimizationfreq( P, omega_ref ), 1e5, 1e7 );

% Fixed-Fixed. 
P2 = fminbnd(@(P) optimizationfreqfixed( P, omega_ref ), 1e5, 1e7 );

% Fixed-Fixed (string)
P3 = fminbnd(@(P) optimizationfreqfixedstring( P, omega_ref ), 1e5, 1e7 );
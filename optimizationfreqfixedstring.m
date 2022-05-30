function [R3] = optimizationfreqfixedstring(P,omega_ref)

n1 = 1;
d1 = 98e-3;           %Diameter of cable [m]
r1 = d1/2;             %Radius of cable [m]
A1 = pi*r1^2;          %Area of cable [m^2]
E1 = 160e9;           %Modulus of elasticity [Pa]
L1 = 177.168;         %Length of cable [m]
rho1 = 51.6;          %Density [kg/m]
I1 = (pi/64)*d1^4;     %Area moment of inertia [m^4]



c = sqrt(P/rho1);

wn = (n1*c*pi)/L1;

R3 = (wn - omega_ref).^2; %Objective function

end 
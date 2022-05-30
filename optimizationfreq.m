function [R1] = optimizationfreq(P,omega_ref)

n1 = 1;
d1 = 98e-3;           %Diameter of cable [m]
r1 = d1/2;             %Radius of cable [m]
A1 = pi*r1^2;          %Area of cable [m^2]
E1 = 160e9;           %Modulus of elasticity [Pa]
L1 = 177.168;         %Length of cable [m]
rho1 = 7830;          %Density [kg/m^3]
I1 = (pi/64)*d1^4;     %Area moment of inertia [m^4]

wn = ((pi^2/L1^2)*sqrt(E1*I1/(rho1*A1))*(n1^4+((n1^2*P*L1^2)/(pi^2*E1*I1))).^(1/2)); 

R1 = (wn - omega_ref).^2; % Objective function 

end 
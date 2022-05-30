function omega = NatFreqFixedFixed(P)
%Created by Luis David Avedaño-Valencia
% Structural parameters
d = 98e-3;           %Diameter of cable [m]
A = pi*(d/2)^2;      %Cross-section area of cable [m^2]
E = 160e9;           %Modulus of elasticity [Pa]
L = 177.168;         %Length of cable [m]
rho = 7830;          %Density [kg/m^3]
I = (pi/64)*d^4;     %Area moment of inertia [m^4]
% P = 2.18e6;

% Set up the boundary conditions
syms a(phi) b(phi) phi
syms f(phi)

lambda = ( P*L^2/(E*I) );
a(phi) = sqrt(  lambda/2 + sqrt( lambda^2/4 + phi ) );
b(phi) = sqrt( -lambda/2 + sqrt( lambda^2/4 + phi ) );

D = [1 0 1 0;
    cosh(a) sinh(a) cos(b) sin(b);
    0 a 0 b;
    a*sinh(a) a*cosh(a) -b*sin(b) b*cos(b)];

% Calculate the eigenvalue equation and its derivatives
f(phi) = simplify( det(D) );                    % Frequency equation 
dfdphi(phi) = simplify( diff( f, phi ) );
r = matlabFunction( f/dfdphi );                 % Ratio of function to derivative (for Newton-Raphson)
f = matlabFunction( f );                        % Eigenvalue function
g = matlabFunction( dfdphi );                   % Eigenvalue function derivative w.r.t. phi

%% Find eigenvalues with the Newton-Raphson method

% Settings of the Newton-Raphson method
tol = 1e-8;
MaxIter = 100;

% Step 1 : Leave the basin of attraction of the zero eigenvalue
x0 = 1;
x0 = LineSearch(f,x0,MaxIter,1e-4,false);

% Initialization of the Newton-Raphson algorithm
x = zeros(1,MaxIter+1);
x(1) = 1.2*x0;
iter = 1;
err = inf;

% Step 2 : Computation loop - Newton-Raphson
while err > tol && iter <= MaxIter
    x(iter+1) = x(iter) - 1e-1*( r(x(iter)) );
    iter = iter + 1;
    err = abs( ( x(iter) - x(iter-1) )/x(iter-1) );
end
x = x(1:iter);
omega = sqrt( x(end)*E*I/(rho*A*L^4) );

end
clc
clear
close all
% Structural parameters
d = 98e-3;           %Diameter of cable [m]
A = pi*(d/2)^2;      %Cross-section area of cable [m^2]
E = 160e9;           %Modulus of elasticity [Pa]
L = 177.168;         %Length of cable [m]
rho = 7830;          %Density [kg/m^3]
I = (pi/64)*d^4;     %Area moment of inertia [m^4]
P = 1030000:1000:1190000; % Axial load P [N]

for i = 1:160

% Set up the boundary conditions
syms a(phi) b(phi) phi
syms f(phi)

lambda = ( P(i)*L^2/(E*I) );
a(phi) = sqrt(  lambda/2 + sqrt( lambda.^2/4 + phi ) );
b(phi) = sqrt( -lambda/2 + sqrt( lambda.^2/4 + phi ) );

D = [1 0 1 0;
    cosh(a) sinh(a) cos(b) sin(b);
    0 a 0 b;
    a.*sinh(a) a.*cosh(a) -b.*sin(b) b.*cos(b)];

% Calculate the eigenvalue equation and its derivatives
f(phi) = simplify( det(D) );                    % Frequency equation 
dfdphi(phi) = simplify( diff( f, phi ) );
r = matlabFunction( f/dfdphi );                 % Ratio of function to derivative (for Newton-Raphson)
f = matlabFunction( f );                        % Eigenvalue function
g = matlabFunction( dfdphi );                   % Eigenvalue function derivative w.r.t. phi

%% Find eigenvalues with the Newton-Raphson method
% Fixed-Fixed beam model 
% Settings of the Newton-Raphson method
tol = 1e-8;
MaxIter = 100;

% Step 1 : Leave the basin of attraction of the zero eigenvalue
x0 = 0.5;
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
omega(i) = sqrt( x(end)*E*I/(rho*A*L^4) );

end

%% Finding eigenvalues for Pinned-Pinned Beam model
n = 1; %First natural frequency

for i=1:160

wn(i) = ((pi^2/L^2)*sqrt(E*I/(rho*A))*(n^4+((n^2*P(i)*L^2)/(pi^2*E*I))).^(1/2));

end 

%% Finding eigenvalues for String model 
rho1 = 51.6; % density in mass per unit length
n = 1;
for i = 1:160
c = sqrt(P(i)/rho1);

wn_string(i) = (n*c*pi)/L;
end 
%% Plotting the 3 physical models 
figure
plot(P(1:160),(omega/(2*pi)))
grid on 
hold on     
plot(P(1:160),(wn/(2*pi)))
grid on 
hold on 
plot(P(1:160),wn_string/(2*pi))
xlabel('Axial Load [N]','Interpreter','latex','FontSize',11)
ylabel('Natural frequency [Hz]','Interpreter','latex','FontSize',11)
legend('Fixed-Fixed','Pinned-Pinned','Fixed-String')

%% Finding the slope and Intercept 
format long 

X = P(1:160)'; % Define X
x = [ones(size(X)), X];

% Fixed-Fixed beam model 
y = omega/(2*pi)';
b1 = regress(y',x)

% Pinned-Pinned Beam model
y1 = wn/(2*pi)';
b2 = regress(y1',x)

% String model 
y2 = wn_string/(2*pi)';
b3 = regress(y2',x)

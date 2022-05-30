%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                      Model based on Physics                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%% Plot of the characteristic equations for the eigenvalues  %%%

clc; clear; close all

d = 98e-3;           %Diameter of cable [m]
r = d/2;             %Radius of cable [m]
A = pi*r^2;          %Area of cable [m^2]
E = 160e9;           %Modulus of elasticity [Pa]
L = 177.168;         %Length of cable [m]
rho = 7830;          %Density [kg/m^3]
rho1 = 51.6;         %Transverse vibrations => density in mass pr. unit length 
I = (pi/64)*d^4;     %Area moment of inertia [m^4]
%P = 10e6;           %Initial guess for axial load
%Axial loads obtained from model updating
%P = 1.0366e6;       %String
%P = 1.1862e6;       %pinned-pinned
P = 1.1652e6;        %Fixed-fixed

%Plot all models together
omega = 0:0.1:100;
[f1,f2,f3,s1,s2] = FixedFixedRoots(omega, E, L, A, I, rho, rho1, P);
figure()
plot(omega/(2*pi),log10(abs(f1)))
hold on
plot(omega/(2*pi),log10(abs(f2)))
hold on 
plot(omega/(2*pi),log10(abs(f3)))
grid on
legend('Fixed-Fixed','Pinned-Pinned','Fixed-string','Location', 'Best')
xlabel('Frequency [Hz]','Interpreter','latex','FontSize',11)
xlim([0 6])

%Plot fixed and pinned model together
figure()
plot(omega/(2*pi),log10(abs(f1)))
hold on
plot(omega/(2*pi),log10(abs(f2)))
grid on
legend('Fixed-Fixed','Pinned-Pinned','Location', 'Best')
title('Fixed-Fixed and Pinned-Pinned','Interpreter','latex','FontSize',13)
xlim([0 6])
xlabel('Frequency [Hz]','Interpreter','latex','FontSize',11)

%Plot each one of the models
figure()
plot(omega/(2*pi),log10(abs(f1)),'Color','[0.8500 0.3250 0.0980]')
grid on
title('Fixed-Fixed','Interpreter','latex','FontSize',13)
xlim([0 2])
xlabel('Frequency [Hz]','Interpreter','latex','FontSize',11)


figure()
plot(omega/(2*pi),log10(abs(f2)),'Color','[0.4940 0.1840 0.5560]')
grid on
title('Pinned-Pinned','Interpreter','latex','FontSize',13)
xlim([0 2])
xlabel('Frequency [Hz]','Interpreter','latex','FontSize',11)

figure()
plot(omega/(2*pi),log10(abs(f3)))
xlim([0 2])
xlabel('Frequency [Hz]','Interpreter','latex','FontSize',11)
title('Fixed-String','Interpreter','latex','FontSize',13)
grid on 


%% %%% Mode shapes plotted for fixed-fixed, pinned-pinned and fixed-string %%%

close all
clc 
clear

%%%%%%% fixed-fixed %%%%%%%

Ndof = 1000;                 %Number of points on the plot
L = 177;                    % Length of beam
E = 160e9;                  % Modulus of elasticity
rho = 7830;                 % Density 
d = 98e-3;                  % Diameter
A = d^2/4 * pi;             % Cross-section area
I = pi/64 * d^4;            % Moment of inertia
n = 4;                      % Number of modes

x = linspace(0,L,Ndof+1);

P = 1.1652e6;               %The axial load used for the three models to plot mode shapes


% Freqency vector for calculation
omega = linspace(0.1,200,5e3);

s1 = sqrt(P/(2*E*I)+sqrt((P/(2*E*I))^2+(rho*A*omega.^2)/(E*I)));             %Root 1 (eq. 8.126 in Rao's book)
s2 = sqrt(-P/(2*E*I)+sqrt((P/(2*E*I))^2+(rho*A*omega.^2)/(E*I)));            %Root 2 (eq. 8.126 in Rao's book)

% Calculating frequency vector
f1 = 2*s1.*s2.*(1-cosh(s1.*L).*cos(s2.*L))+(s1.^2-s2.^2).*sinh(s1.*L).*sin(s2.*L);

% Natural frequencies are found as local minima
indx = islocalmin(log10(abs(f1)));

% Extracting the analytical natural frequencies
Freq = omega(indx) / (2*pi);

% Calculating modeshapes
Psi = zeros(Ndof+1,n); 
s1 = s1(indx);
s2 = s2(indx);
for i = 1:n
    % Constants used for shape function
    C1 = 1;
    C3 = -C1;
    C4 = C1 * (cos(s2(i)*L)/sin(s2(i)*L));
    C2 = (-(C1*cosh(s1(i)*L)-C1*cos(s2(i)*L)+C4*sin(s2(i)*L)))/(sinh(s1(i)*L));
    
    % Calculating modeshapes for each mode
    Psi(:,i) = C1*cosh(s1(i)*x)+C2*sinh(s1(i)*x)+C3*cos(s2(i)*x)+C4*sin(s2(i)*x);
end

j = Psi(:,1)/max(Psi(:,1));
k = Psi(:,2)/max(Psi(:,2));
l = Psi(:,3)/max(Psi(:,3));
m = Psi(:,4)/max(Psi(:,4));

Psi = [j, k, l, m];

% Plotting modeshapes 
figure
for i = 1:n
    subplot(2,2,i)
    plot(x, Psi(:,i) ,'LineWidth',2)


    xlabel('Length [m]','Interpreter','latex','FontSize',11)
    ylabel('Amplitude','Interpreter','latex','FontSize',11)
    title(['Mode ',num2str(i)],'Interpreter','latex','FontSize',13)
    sgtitle('Fixed-Fixed','Interpreter','latex','FontSize',13) 
    grid on
    xlim([0 180])   
    hold on
end


%%%%%%% pinned-pinned %%%%%%%

%%%%%            Natural frequencies and mode shapes             %%%%%%%
% Freqency vector for calculation
omega = linspace(0.1,200,5e3);

s1 = sqrt(P/(2*E*I)+sqrt((P/(2*E*I))^2+(rho*A*omega.^2)/(E*I)));             %Root 1 (eq. 8.126 in Rao's book)
s2 = sqrt(-P/(2*E*I)+sqrt((P/(2*E*I))^2+(rho*A*omega.^2)/(E*I)));            %Root 2 (eq. 8.126 in Rao's book)

% Calculating frequency vector
f1 = sinh(s1.*L).*sin(s2.*L); %Pinned-Pinned (eq. E6 in Rao's book)

% Natural frequencies are found as local minima
indx = islocalmin(log10(abs(f1)));

% Extracting the analytical natural frequencies
Freq = omega(indx) / (2*pi);

% Calculating modeshapes
Psi = zeros(Ndof+1,n);  
s1 = s1(indx);
s2 = s2(indx);
for i = 1:n
    % Constants used for shape function
    C1 = 1;
    C3 = -C1;
    C4 = C1 * (cos(s2(i)*L)/sin(s2(i)*L));
    C2 = (-(C1*cosh(s1(i)*L)-C1*cos(s2(i)*L)+C4*sin(s2(i)*L)))/(sinh(s1(i)*L));

    
    % Calculating modeshapes for each mode
    Psi(:,i) = -(C1*cosh(s1(i)*x)+C2*sinh(s1(i)*x)+C3*cos(s2(i)*x)+C4*sin(s2(i)*x));
end

j = Psi(:,1)/max(Psi(:,1));
k = Psi(:,2)/max(Psi(:,2));
l = -Psi(:,3)/max(Psi(:,3));
m = Psi(:,4)/max(Psi(:,4));

Psi = [j, k, l, m];


% Plotting modeshapes 

for i = 1:n
    subplot(2,2,i)
    plot(x, Psi(:,i) ,'LineWidth',2)

    xlabel('Length [m]','Interpreter','latex','FontSize',11)
    ylabel('Amplitude','Interpreter','latex','FontSize',11)
    title(['Mode ',num2str(i)],'Interpreter','latex','FontSize',13)
    sgtitle('Pinned-Pinned','Interpreter','latex','FontSize',13) 
    grid on
    xlim([0 180])
    hold on
end


%%%%%%% fixed-string %%%%%%%

E = 160e9;                  % Modulus of elasticity
rho = 51.6;                 % Transverse vibrations => density in mass pr. unit length 
d = 98e-3;                  % Diameter
A = d^2/4 * pi;             % Cross-section area
I = pi/64 * d^4;            % Moment of inertia
n = 1:4;                    % Number of modes


c = (P/rho)^(1/2);


Psi = zeros(Ndof+1,4);  
for i = 1:4

    omega(i) = ((n(i)*c*pi)/L);

    Psi(:,i) = sin((omega(i).*x)./c);

end

% Plotting modeshapes 
for i = 1:4
    subplot(2,2,i)
    %plot(((1:1:Ndof+1)-1), Psi(:,i) ,'LineWidth',2)
    plot(x, Psi(:,i) ,'LineWidth',2)

    xlabel('Length [m]','Interpreter','latex','FontSize',11)
    ylabel('Amplitude','Interpreter','latex','FontSize',11)
    title(['Mode ',num2str(i)],'Interpreter','latex','FontSize',11)
    sgtitle('Mode shapes','Interpreter','latex','FontSize',13) 
    grid on
    legend('Fixed-Fixed','Pinned-Pinned','Fixed-String','Interpreter','latex','FontSize',10)
    xlim([0 180])
    
end

function [f1,f2,f3,s1,s2] = FixedFixedRoots(omega, E, L, A, I, rho, rho1, P)

s1 = sqrt(P/(2*E*I)+sqrt((P/(2*E*I))^2+(rho*A*omega.^2)/(E*I)));    %Root 1 (eq. 8.126 in Rao's book)

s2 = imag(sqrt(P/(2*E*I)-sqrt((P/(2*E*I))^2+(rho*A*omega.^2)/(E*I))));    %Root 2 (eq. 8.126 in Rao's book)

f1 = 2*s1.*s2.*(1-cosh(s1.*L).*cos(s2.*L))+(s1.^2-s2.^2).*sinh(s1.*L).*sin(s2.*L);   %Fixed-Fixed (eq. 6.16 in the book "Structural Stability and Vibration")

f2 = sinh(s1.*L).*sin(s2.*L); %Pinned-Pinned (eq. E6 in Rao's book)

c = (P/rho1)^(1/2);
f3 = sin((omega*L)/c);  %Fixed-String

end



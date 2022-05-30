function [R2] = optimizationfreqfixed(P,omega_ref)

omega = NatFreqFixedFixed(P);
R2 = (omega - omega_ref).^2; %Objective function

end 
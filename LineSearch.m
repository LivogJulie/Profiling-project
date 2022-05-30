function [x,n] = LineSearch(f,delta,kmax,tol,PlotFlag)
%--------------------------------------------------------------------------
% This function performs phases 0, 1, 2 and 3 of the line search algorithm
% using the Golden section algorithm
% Inputs :      f       -> Function handle
%               delta   -> Initial step size
%               kmax    -> Maximum number of iterations
%               tol     -> Termination tolerance for the Golden Section method
% Outputs :     x       -> Optimal point (minimum of f(x))
%               n       -> Total number of cost function evaluations
%
% Created by : Andriana Georgantopoulou - December 2021
%--------------------------------------------------------------------------

if nargin < 5
    PlotFlag = false;
end

% Defining computation variables
gamma = ( 1 + sqrt(5) ) / 2;                                                % Golden ratio
tau = gamma - 1;                                                            % Golden ratio - 1
n = 0;                                                                      % Number of function evaluations

% Phase 0: Adjustment of initial point
k = 0;
f0 = f(0);
n = n+1;
while f0 <= f(delta) && k<=kmax
    delta = delta/10;           
    k=k+1; n=n+1;
end

% Phase 1: Golden Section Phase 1 (bracketing the minimum)
k = 1;
alpha = delta;
ff = [f0 f(alpha)];
while ff(k+1)<ff(k) && k<=kmax

    alpha(k+1) = alpha(k) + gamma.^(k+1)*delta;                        % New alpha
    ff(k+2) = f(alpha(k+1));                                % Function value at new alpha
    
    n = n+1;                                                % Increase function evaluation count
    k = k+1;                                                % Increase iteration count

end
alpha_u = alpha(k);
if k>2
    alpha_l = alpha(k-2);
else 
    alpha_l = 0;
end

% Phase 2: Closing the interval
I = alpha_u - alpha_l;                                                      % Initial interval size
k = 0;                                                                      % Number of iterations

% Main iteration loop - Finishes if the interval size is smaller than the
% tolerance  or if the number of iterations reaches the limit
while I >= tol && k <= kmax                                                  
    
    % Get new interval candidate values
    alpha_a = alpha_l + (1-tau)*I;                                          % New lower limit candidate
    alpha_b = alpha_l + tau*I;                                              % New higher limit candidate
    f_a = f(alpha_a);                                                       % Function value in alpha_a
    f_b = f(alpha_b);                                                       % Function value in alpha_b
    
    n = n+2;                                                                % Increase number of function evaluations
    k = k+1;                                                                % Increase iteration number

    % Select new interval limits
    if f_a < f_b
        alpha_u = alpha_b;
    elseif f_a > f_b
        alpha_l = alpha_a;
    elseif f_a == f_b
        alpha_l = alpha_a;
        alpha_u = alpha_b;
    end

    I = alpha_u - alpha_l;                                                  % New interval size

    if PlotFlag                                                             % Plot current iteration progress
        semilogy(k,I,'ob')
        hold on
    end

end

if PlotFlag
    plot([0 k],tol*[1 1],'--k')
    xlabel('No. Iterations - Phase 2')
    ylabel('Interval size')
    grid on
end

% Phase 3: Quadratic interpolation
alpha_i = alpha_l + 0.5*I;

fU = f(alpha_u);
fL = f(alpha_l);
fI = f(alpha_i);
n=n+1;

a2 = 1/(alpha_u-alpha_i)*((fU-fL)/(alpha_u-alpha_l) - (fI-fL)/(alpha_i-alpha_l));
a1 = ((fI-fL)/(alpha_i-alpha_l)) - a2*(alpha_l+alpha_i);

if a2<=0 || a1==0
    x = alpha_i;                                                            % Solution before quadratic interpolation
else
    x=-a1/(2*a2);
end

% alpha = [alpha_l alpha_l+0.5*I alpha_u]';                                   % Interpolation points
% M = [ones(3,1) alpha alpha.^2];                                             % Quadratic function values
% F = zeros(3,1);
% for i=1:3
%     F(i) = f(alpha(i));                                                     % Evaluate the function on the interpolation points
% end
% n = n+3;
% 
% a = pinv(M)*F;                                                              % OLS with pseudo-inverse
% if a(2) > 0
%     x = -a(2)/(2*a(3));                                                     % Minimum of the obtained quadratic function
% end

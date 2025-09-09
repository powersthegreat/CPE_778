clc; clear;

% Add path for custom functions (if needed)
addpath('../../functions');

% Define function and derivatives as function handles
f  = @(x) (x - 2).^2;        % Objective function
F  = @(x) 2*(x - 2);         % First derivative
H  = @(x) 2;                 % Second derivative

% Initial guess
x0 = 10;

% Run NaiveNewton
[x, fun, dfun, xk] = naiveNewton(f, F, H, x0);

% Display results
fprintf('\nFinal solution:\n');
fprintf('x   = %.8f\n', x);
fprintf('f(x) = %.8f\n', fun);
fprintf('f''(x) = %.8f\n', dfun);

% Plot convergence history
figure;
plot(xk, 'o-','LineWidth',1.5);
xlabel('Iteration');
ylabel('x value');
title('Convergence of Newton''s Method');
grid on;
